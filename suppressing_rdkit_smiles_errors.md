# Suppressing RDKit SMILES Parse Errors

A copy-paste fix for any repo that parses untrusted SMILES with RDKit and wants
to silence the noisy `[HH:MM:SS] SMILES Parse Error: ...` lines RDKit prints to
stderr when `Chem.MolFromSmiles()` gets garbage.

## The fix (one line)

```python
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
```

That's it. After this call, `Chem.MolFromSmiles('this is not smiles')` returns
`None` quietly instead of spamming stderr with a multi-line parse trace.

### Where to put it

Put it **as early as possible**, before the first `Chem.MolFromSmiles` /
`Chem.MolFromSmarts` / `rdMolDescriptors` call in your program — ideally right
after the rdkit import at the top of your entrypoint (or your CLI's `main()`),
and before any module that parses SMILES gets exercised.

```python
# main.py  (or wherever your program starts)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   # silence invalid-SMILES stderr noise

import sys
sys.path.insert(0, 'code')      # then import your helpers
from MolPropOp import grow_cycle, replace_groups
```

It is global and process-wide: rdkit uses a single logger, so disabling it once
affects every subsequent call in the process, including calls made inside
imported libraries. You do **not** need to repeat it in every file.

## Why this is the right level of suppression

`Chem.MolFromSmiles` is designed to *not* raise on bad input — it returns `None`
and logs the parse error. The return-`None` contract is what you actually want
for untrusted input (you branch on `if mol is None`). The stderr logging is the
only annoying part, and `rdApp.*` is exactly the logger that owns it. Disabling
it keeps the safe `None`-returning behaviour and just drops the chatter.

You are **not** hiding a crash or a real exception — `MolFromSmiles` never raised
to begin with. You're only turning off a diagnostic print.

## Narrower alternative

If you'd rather keep other RDKit application logs and only mute the SMILES
parser, target the specific logger instead of the whole `rdApp` tree:

```python
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.smiles')
```

`'rdApp.*'` is the blunt instrument (mutes everything under the application
logger, including sanitization warnings); `'rdApp.smiles'` is the scalpel. In
practice the SMILES parser is the only loud one during batch parsing, so either
works — use `'rdApp.smiles'` if you want to still see other RDKit warnings.

## Re-enabling for debugging

When something genuinely won't parse and you want the trace back (e.g. a CLI
debug flag), flip it on again:

```python
RDLogger.EnableLog('rdApp.*')   # or 'rdApp.smiles'
```

A typical CLI pattern:

```python
parser.add_argument('--rdkit-verbose', action='store_true',
                    help='Re-enable RDKit stderr (SMILES Parse Error) logs. '
                         'They are silenced by default since invalid input is '
                         'surfaced to the caller as a None return / a skipped row.')
# ...
if args.rdkit_verbose:
    RDLogger.EnableLog('rdApp.*')
```

## Pair it with surfacing the failure to the caller

Silencing the log is only half the fix. If your code loops over candidate
SMILES, make sure a `None` result is *reported back* (logged at your level, or
appended as a skipped entry) instead of silently dropped — otherwise the consumer
can't tell "0 results because all invalid" from "0 results because nothing
matched". Example:

```python
mol = Chem.MolFromSmiles(smi)
if mol is None:
    results.append((smi, 'invalid SMILES, skipped'))
    continue
```

Without this, silencing the stderr noise can hide the fact that every input was
garbage. Silence the logger, but never swallow the signal.

## Gotchas

- **Call it before the first parse.** If a module parses SMILES at import time,
  the disable must run before that module is imported.
- **It does not suppress Python tracebacks.** If you're calling something that
  *raises* on bad SMILES (e.g. `Chem.MolFromSmiles(smi, sanitize=True)` with a
  custom error handler, or `Chem.MolFromSmiles` inside a block that re-raises),
  wrap that in your own `try/except`. `RDLogger.DisableLog` only touches RDKit's
  logger, not Python exceptions.
- **Other RDKit APIs that log**: `MolFromSmarts`, `MolFromMolBlock`,
  `rdMolDescriptors.*`, and the reaction parser all log under `rdApp.*` too —
  disabling the tree mutes all of them, which is usually what you want for batch
  jobs.
- **Subprocess/pipe note**: the parse errors go to the process's *stderr*. If
  you've redirected `2>&1` into a log file (common in CLI runs), that's where
  they accumulate — this fix removes them at the source instead.

## One-file self-test

Drop this into a scratch file to confirm the behaviour before wiring it into
your repo:

```python
from rdkit import RDLogger, Chem

print('--- suppressed ---')
RDLogger.DisableLog('rdApp.*')
print('returned:', Chem.MolFromSmiles('this is not smiles'))   # None, no stderr

print('--- re-enabled ---')
RDLogger.EnableLog('rdApp.*')
print('returned:', Chem.MolFromSmiles('still not smiles'))     # None + stderr trace
```

Expected: the suppressed block prints only `returned: None`; the re-enabled
block also prints the `SMILES Parse Error` lines.