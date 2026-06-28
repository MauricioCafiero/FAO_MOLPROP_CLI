# FAO_MOLPROP_CLI

Adversarial molecule optimization on the command line. An **Ollama** main model
designs molecules to bind a protein target, calling chemistry + docking tools to
propose and score them; an **adversary** model (OpenAI or Anthropic) critiques
each proposal; the main model refines until it signals `Done`. Every step is
written to a timestamped Markdown report under `results/`, with a JSON messages
sidecar written alongside it so a stopped run can be resumed.

`molopt.py` is a single-file CLI port of `Ollama_MolOpt.ipynb`, which is kept in
the repo as the reference/source of truth for the prompts.

## What it does

Headless loop, one target at a time:

1. The main (Ollama) model is seeded with a starting molecule/score list
   (`code/adversarial_set.md`) and reasoning instructions.
2. It reasons about the data and calls **chemistry tools** (grow a scaffold,
   swap groups, dock, check drug-likeness, …) to explore chemical space.
3. Each proposal is sent to the **adversary**, which finds flaws in the
   reasoning / estimated docking scores and suggests modifications.
4. The main model revises in light of the critique. The loop ends when the main
   model replies `Done` (robustly detected: a trailing standalone `Done` after a
   final proposal counts) or when `--max-turns` is hit.

Docking is done with **dockstring** against a receptor (default target HMGCR);
the receptors `HMGCR_dude_receptor.pdb` / `HMGCR_dude_receptor_2.pdb` live at the
repo root. Docking is the runtime bottleneck, not the LLM.

## Architecture

- **Main model** — Ollama, function-calling. Hosted default `https://ollama.com`
  with a bearer token, or a local daemon via `--ollama-host
  http://localhost:11434`. `--think`/`--no-think` auto-enables thinking for
  models that support it (the built-in set all do; `NO_THINK_MODELS` is empty).
- **Adversary** — `--adversary openai|anthropic`. OpenAI uses the Responses API
  (`client.responses.create`); Anthropic uses the Messages API
  (`client.messages.create`). SDKs are lazy-imported so a missing dep doesn't
  block the other. The Anthropic client passes an explicit
  `base_url="https://api.anthropic.com"` (overridable via `--anthropic-base-url`)
  so it bypasses any ambient `ANTHROPIC_BASE_URL` pointing at a local proxy.
- **Chemistry tools** — reused unchanged from `code/`:
  - `MolPropOp.py`: `grow_cycle`, `replace_groups`, `make_random_list`,
    `related`, `lipinski`
  - `docking_module.py`: `dock_and_get_interacting_residues`,
    `calculate_SAS_and_NP`, plus the shared mutable `scoring_args`
    (`[cpu_count, protein]`) and the NumPy-2 `np.in1d` shim.
- **Robustness** — auto-retry on Ollama connection errors (`--max-retries`);
  tool-call exceptions are fed back to the model for self-correction;
  `--max-tool-calls` caps per-turn tool rounds (once hit, a text-only summary is
  forced — a stuck model can't loop forever); invalid substituents are surfaced
  to the model as `invalid SMILES, skipped` entries so it stops retrying them.

## Repository layout

```
molopt.py                     # the CLI (entrypoint)
test_models.py                # smoke-test candidate Ollama cloud models
Ollama_MolOpt.ipynb           # original notebook (reference, not deleted)
new_models_to_add.txt         # candidate cloud model names (with -cloud suffix)
.env.example                  # copy to .env and fill in API keys
requirements.txt              # full deps (base + chem stack + LLM SDKs)
requirements-base.txt         # minimal build deps (setuptools/wheel/six)
requirements-oddt.txt         # optional: oddt (heavier scientific stack)
suppressing_rdkit_smiles_errors.md  # drop-in fix to silence RDKit parse logs
code/
  adversarial_set.md          # starting molecule/score list fed to the model
  MolPropOp.py                # molecular-property operations (grow/replace/...)
  docking_module.py           # dockstring docking + SAS/NP scoring
HMGCR_dude_receptor_2.pdb      # the receptor actually used by the code (hardcoded in docking_module.py)
HMGCR_dude_receptor.pdb        # older copy, kept as a backup; not referenced by the code
results/                       # timestamped session reports (.md) + JSON message sidecars (.json)
```

## Install

Install in this **specific order** — it matters because `oddt` has been
unmaintained for years and breaks the rest of the stack if installed normally.
Install the build tools first, then `oddt` *without build isolation* (so it
reuses the setuptools/wheel just installed instead of trying to fetch its own,
ancient build deps), then everything else last.

```bash
python3 -m venv fao-env
source fao-env/bin/activate

# 1) build tools first (oddt's --no-build-isolation install needs them present)
pip install -r requirements-base.txt        # setuptools, wheel, six

# 2) oddt, installed WITHOUT build isolation (unmaintained; would otherwise fail
#    to build and/or drag in incompatible deps). Required by docking_module.py.
pip install --no-build-isolation -r requirements-oddt.txt   # oddt

# 3) all other packages last (rdkit, dockstring, pyscf, openbabel-wheel,
#    ollama, openai, anthropic, ...)
pip install -r requirements.txt
```

> **Why the order matters:** `oddt` is unmaintained and its build is fragile. If
> you `pip install -r requirements.txt` first and then `oddt`, or let `oddt`
> build in isolation, it can pull an incompatible NumPy / fail to compile, which
> breaks `docking_module.py`. Installing build tools → `oddt` with
> `--no-build-isolation` → the rest avoids that. `docking_module.py` also
> applies an `np.in1d -> np.isin` shim at import time for NumPy 2.x, since oddt
> 0.7 otherwise crashes on NumPy 2.
>
> Note: `openbabel-wheel` and `rdkit` are pip-installable wheels on macOS/manylinux.

## Configure

Copy `.env.example` to `.env` and fill in keys, **or** export them in your shell
(`source ~/.zshrc` if they live there). Real environment variables always override
`.env`, which is only read from the directory you run `molopt.py` from.

| Variable          | Used by                          | Required?                          |
|-------------------|----------------------------------|------------------------------------|
| `OLLAMA_API_KEY`  | main model (hosted endpoint)     | yes for `https://ollama.com`; no for local daemon |
| `OPENAI_API_KEY`  | `--adversary openai`             | only if using the OpenAI adversary |
| `ANTHROPIC_API_KEY` | `--adversary anthropic`        | only if using the Anthropic adversary |

## Usage

```bash
# list the built-in Ollama model set, then run with the default model
python3 molopt.py --list-models
python3 molopt.py --protein HMGCR

# specific cloud model (API name, no '-cloud' suffix) with think forced on
python3 molopt.py --model gemma4:31b --think

# Anthropic adversary instead of the default OpenAI one
python3 molopt.py --model qwen3.5:397b \
  --adversary anthropic --adversary-model claude-haiku-4-5-20251001

# quick chemistry-stack check with no LLM keys
python3 molopt.py --self-test
```

See `python3 molopt.py --help` for the full option list and an examples block.

### Key flags

| Flag | Default | Purpose |
|------|---------|---------|
| `--protein` | `HMGCR` | Docking target. |
| `--model` | `deepseek-v4-pro` | Ollama main model (API name, no `-cloud` suffix). |
| `--adversary` | `openai` | `openai` (default `gpt-5.2`) or `anthropic` (default `claude-haiku-4-5-20251001`). |
| `--max-turns` | `20` | Safety cap on adversary↔main turns. |
| `--max-tool-calls` | `12` | Per-turn cap on tool-calling rounds; forces a text summary once hit. |
| `--ollama-host` | `https://ollama.com` | Use `http://localhost:11434` for a local daemon. |
| `--context-file` | `code/adversarial_set.md` | Starting molecule/score list. Ignored when `--resume` is given. |
| `--resume` | off | Resume a prior run from its JSON sidecar; skips context seeding + the initial model turn and continues from the last assistant message. |
| `--results-dir` | `./results` | Where timestamped reports + sidecars are written. |
| `--trace` | off | Write the full thinking + tool-call trace into the results `.md`. |
| `--rdkit-verbose` | off | Re-enable RDKit `SMILES Parse Error` stderr (silenced by default). |
| `--quiet` | off | Suppress thinking/content/tool prints to stdout. |
| `--self-test` | — | Run the chemistry tools directly with no LLM keys, then exit. |
| `--list-models` | — | Print the built-in model list and exit. |

## Models

The built-in Ollama set (see `--list-models`) is the cloud catalog with the
`-cloud`/`:cloud` suffix dropped — that suffix is a catalog label; the **chat API
name drops it** (e.g. `gemma4:31b-cloud` → `gemma4:31b`, `glm-5.2:cloud` →
`glm-5.2`). `new_models_to_add.txt` lists the catalog names; `molopt.py` uses the
stripped API names. All five current models support thinking + tool-calling
simultaneously (verified by `test_models.py`), so `NO_THINK_MODELS` is empty; the
old model list and old no-think set are kept as comments for backwards
compatibility.

## Outputs

Each run writes `results/<model>_<protein>_<timestamp>.md` containing the turn
headers, the main model's proposals, and the adversary feedback. By default only
terse section headers + text go to the `.md`; the full thinking/tool-call trace
goes to stdout (use `--trace` to include it in the `.md` too). The final line
reports the terminal status: `# Session end: Done` (model signalled Done) or a
`MAX_TURNS` message.

### Example: GLM-5.2 + Claude adversary optimizing HMGCR

A full adversarial run (the command used for the example below):

```bash
python3 -u molopt.py --protein HMGCR --model glm-5.2 \
  --adversary anthropic --adversary-model claude-haiku-4-5-20251001 \
  --max-turns 2 --max-tool-calls 3
```

`--max-tool-calls 3` is the practical ceiling: each `grow_cycle` docks ~45
molecules, so the default `12` would mean ~540 dockings per turn (~75 min/turn).
`3` keeps a turn to ~135 dockings (~15-20 min) while still letting the model
explore before the adversary critiques. (Once the cap is hit the model is
nudged to emit a text summary — *without* saying `Done` — so the adversary
loop still runs.)

The run completed two full adversary↔main cycles (initial model turn, then
Turn 1 and Turn 2), each adversary critique driving a refinement:

| Phase | What happened |
|-------|---------------|
| Initial turn | Model explores the coumarin-flavone scaffold via `grow_cycle`; reaches ~-9.3 with bis-phenyl carboxylates (LogP ~5.6, **fails Lipinski**) |
| Adversary 1 | "Analysis of Proposed Molecules" — flags the over-lipophilic phenyl-rich hits, suggests heteroaromatic replacements |
| Turn 1 | Model applies the suggestion: replaces pendant phenyl with 3-pyridinyl |
| Adversary 2 | "Critical Analysis of Proposed Molecules" — confirms the pyridine gain, pushes for F substitution + validation |
| Turn 2 | Model lands the breakthrough: **-9.5** with LogP **4.44** (now **passes Lipinski**) |

Best molecules (dockstring units, more negative is better; Mac scores — see
the macOS caveat below):

| # | SMILES | Score | QED | LogP | Lipinski | SAS |
|---|-------|:---:|:---:|:---:|:---:|:---:|
| 1 (di-pyridine + F) | `O=c1cc(-c2c(c7cccnc7)cc(F)cc2)oc2c(c7cccnc7)ccc(C(C(=O)[O-]))c12` | **-9.5** | ~0.40 | ~4.0 | ✅ | pending |
| 2 (di-pyridine + Me) | `O=c1cc(-c2c(c7cccnc7)cc(C(C))cc2)oc2c(c7cccnc7)ccc(C(C(=O)[O-]))c12` | **-9.5** | 0.371 | 4.44 | ✅ | 2.98 |
| 3 (phenyl + Me) | `O=c1cc(-c2c(c7ccccc7)cc(C(C))cc2)oc2cccc(C(C(=O)[O-]))c12` | **-9.0** | 0.523 | 3.98 | ✅ | 2.66 |

The SAR the loop converged on: keep the chromone carboxylate anchor
(`C(C(=O)[O-])`), and **replace pendant phenyls with 3-pyridinyl** — the
nitrogen adds H-bond acceptors and polar surface area, improving *both* binding
(-9.3 → -9.5) and drug-likeness (LogP 5.65 → 4.44). Position matters:
3-pyridinyl (-9.5) >> 4-pyridinyl (-8.9) >> 2-pyridinyl (-8.3).

Results for this run: `results/glm-5.2_HMGCR_2026-06-28_12-44-57.md`
(+ `.json` messages sidecar — see [Resuming a run](#resuming-a-run)).

### Resuming a run

Every run also writes a JSON messages sidecar next to the `.md` report (same
timestamp stem, `.json`), updated incrementally after each turn (atomically, so
a killed run still leaves a usable sidecar of the last completed turn). Use
`--resume <sidecar.json>` to pick up a stopped run from where it left off
instead of re-seeding from `code/adversarial_set.md`:

```bash
# resume the example run from its sidecar, for one more refinement turn
python3 -u molopt.py --resume results/glm-5.2_HMGCR_2026-06-28_12-44-57.json \
  --protein HMGCR --model glm-5.2 --adversary anthropic \
  --adversary-model claude-haiku-4-5-20251001 --max-turns 1
```

On resume the conversation is seeded from the sidecar (the model already knows
what it proposed and scored), context-file seeding and the initial model turn
are skipped, and refinement continues from the last assistant message. A fresh
timestamped `.md` + sidecar are written for the resumed session; the source
sidecar is never modified.

## Smoke-testing models

`test_models.py` is a fast standalone check (no chem-stack import) that, per
model: resolves the API name (stripped first, raw fallback), tests tool-use,
tests `think=True` + tools together, and runs one short adversary turn. Keys are
resolved the same way as `molopt.py` (CLI flags → env → `.env`).

```bash
source ~/.zshrc   # for the API keys
python3 -u test_models.py                       # all models in new_models_to_add.txt
python3 -u test_models.py --only glm-5.2:cloud  # just one
```

## Notes & caveats

- **macOS docking**: dockstring prints a `DockstringWarning` that Mac scores may
  not match the Linux baselines in the DOCKSTRING paper. Scores are usable for
  relative ranking within a run; don't compare them directly to paper baselines.
- **RDKit stderr**: invalid-SMILES `SMILES Parse Error` lines are silenced by
  default (`RDLogger.DisableLog('rdApp.*')`); failures are surfaced to the model
  as `invalid SMILES, skipped` instead. See
  [`suppressing_rdkit_smiles_errors.md`](suppressing_rdkit_smiles_errors.md) for a
  drop-in copy of that fix for other repos. `--rdkit-verbose` re-enables the log.
- **`ANTHROPIC_BASE_URL` gotcha**: the Anthropic adversary routes to
  `api.anthropic.com` explicitly so a locally-set `ANTHROPIC_BASE_URL` proxy
  (which 404s on `claude-*` names) is bypassed. Use `--anthropic-base-url` to
  point at your own proxy if you actually want one.
- **Keys per shell**: each fresh shell does not auto-load `.zshrc` exports, so
  `source ~/.zshrc` (or use `.env`) before launching.

## TODO / later polish

- Suppress the per-molecule dockstring `DockstringWarning` to clean up `run.log`.