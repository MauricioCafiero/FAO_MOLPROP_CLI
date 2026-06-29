#!/usr/bin/env python3
"""
verify_results.py - Verify a molopt.py run's final molecules and log them to sqlite.

Pipeline:
  1. Parse a results .md written by molopt.py.
  2. Take the LAST model-response block (`# Model response:`, falling back to
     `# Initial model response:` if there were no adversary turns) and extract
     every RDKit-valid SMILES it contains.
  3. For each candidate, recompute the five metrics by REUSING the project's
     existing helper functions (not reimplementing them):
       - docking score  <- docking_module.scoring_function
       - QED, aLogP     <- MolPropOp.lipinski   (aLogP = RDKit QED.properties LogP)
       - SAS, NP        <- docking_module.calculate_SAS_and_NP
  4. Insert novel molecules into a sqlite DB. Novelty key = canonical SMILES +
     InChIKey. Molecules already present are skipped (left untouched).

The md header records the protein and main model, so the docking target is set
automatically (it mutates the shared `scoring_args` list exactly the way
molopt.py does at runtime -- see molopt.py:537-538).

Usage:
  python3 verify_results.py results/<run>.md
  python3 verify_results.py                       # auto-pick newest results/*.md
  python3 verify_results.py <run>.md --dry-run    # list extracted SMILES, no docking/DB
  python3 verify_results.py <run>.md --db molecules.sqlite --min-heavy-atoms 5
"""

import os
import re
import sys
import glob
import sqlite3
import argparse
import contextlib
from datetime import datetime


@contextlib.contextmanager
def _quiet_stdout():
    """Suppress stdout while calling the reused helpers.

    lipinski() prints "lipinski tool"/"=====" and npscorer prints "reading NP
    model ..."/"model in" -- useful in the notebook, but noise here that breaks
    up the results table. Stderr (RDKit is already silenced) is left alone.
    """
    saved = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        yield
    finally:
        sys.stdout.close()
        sys.stdout = saved

# --- Path + NumPy-2 shim setup (must run before importing the helpers) ------
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, 'code'))

# ODDT 0.7 calls np.in1d, removed in NumPy 2.x. np.isin is a drop-in. docking_module
# shims this too, but set it first so the import order is irrelevant.
import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin

from rdkit import Chem, RDLogger
# We MolFromSmiles many non-SMILES tokens while extracting; silence the noise.
RDLogger.DisableLog('rdApp.*')

# Existing helper functions -- reused, not reimplemented.
from docking_module import scoring_function, calculate_SAS_and_NP, scoring_args
from MolPropOp import lipinski


# --- md parsing --------------------------------------------------------------

# A top-level section header written by molopt.py. Matched with `^# ` (single
# hash + space) so the model's own `##`/`###` sub-headers inside a response do
# NOT end the block -- only molopt's known section headers do.
_SECTION_START_RE = re.compile(
    r'^# (Initial model response|Model response|Adversary feedback|Session end'
    r'|Resumed from sidecar|Last assistant text at resume|Trace'
    r'|protein|Adversarial Design Session)\b',
    re.M)


def parse_header(md_text):
    """Pull protein and main-model name out of the md header line.

    Header looks like:
      # protein: HMGCR | main model: glm-5.2 (think=True) | adversary: ...
    """
    protein, model = None, None
    for line in md_text.splitlines():
        if line.startswith('# protein:') or line.startswith('#protein:'):
            m = re.search(r'protein:\s*(\S+)', line)
            if m:
                protein = m.group(1)
            m = re.search(r'main model:\s*(\S+)', line)
            if m:
                model = m.group(1)
            break
    return protein, model


def last_model_response_block(md_text):
    """Return the body text of the last model-response section in the md.

    Recognises `# Initial model response:` and `# Model response:`. The body
    runs from the line after the header up to (but not including) the next
    molopt top-level section. Returns '' if the run ended before any model turn
    produced text (e.g. a killed run with only the header written).
    """
    starts = [(m.end(), m.group(1))
              for m in _SECTION_START_RE.finditer(md_text)
              if m.group(1) in ('Initial model response', 'Model response')]
    if not starts:
        return ''
    start = starts[-1][0]
    nxt = _SECTION_START_RE.search(md_text, start)
    end = nxt.start() if nxt else len(md_text)
    return md_text[start:end].strip()


# --- SMILES extraction -------------------------------------------------------

# Candidate-SMILES regex. Backtick and '|' are intentionally NOT in the class,
# so they act as boundaries -- the pattern pulls SMILES straight out of `...`
# code spans and markdown table cells without needing to parse the table.
# False positives (e.g. "-9.50") are filtered downstream by RDKit validation.
_SMILES_PATTERN = re.compile(r'[CHONFClBrISPKacnosp0-9@+\-\[\]\(\)\/.=#$%]{5,}')


def extract_smiles(text, min_heavy_atoms):
    """Return [(original_smiles, canonical_smiles, mol)] for valid, novel-by-form SMILES.

    Pulls candidate tokens with the SMILES regex, then keeps a candidate only if
    RDKit parses it, the mol has >= min_heavy_atoms heavy atoms, and it isn't a
    bare inorganic ion. Dedupes by canonical SMILES.
    """
    out = []
    seen_canonical = set()
    for m in _SMILES_PATTERN.finditer(text):
        cand = m.group(0)
        # Trim a leading/trailing '.' (disconnected-structure separator) which
        # the regex can grab from prose like "...c12. The" without changing the
        # molecule; never strip it from the middle of a valid SMILES.
        cand = cand.strip('.')
        if not cand or len(cand) < 5:
            continue
        mol = Chem.MolFromSmiles(cand)
        if mol is None:
            continue
        hac = mol.GetNumHeavyAtoms()
        if hac < min_heavy_atoms:
            continue
        # Drop bare ions / salts the model didn't really propose (NaCl, K+, ...).
        if not any(a.GetIsAromatic() or a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
            continue
        canon = Chem.MolToSmiles(mol)
        if canon in seen_canonical:
            continue
        seen_canonical.add(canon)
        out.append((cand, canon, mol))
    return out


# --- metric parsing (reuse the helpers' string outputs) ----------------------

def _parse_lipinski(smiles):
    """Call lipinski([smiles]) and pull out QED and aLogP (RDKit QED.properties LogP)."""
    s = lipinski([smiles])
    qed = alogp = None
    m = re.search(r'QED:\s*([-\d.]+)', s)
    if m:
        qed = float(m.group(1))
    m = re.search(r'LogP:\s*([-\d.]+)', s)
    if m:
        alogp = float(m.group(1))
    return qed, alogp


def _parse_sas_np(smiles):
    """Call calculate_SAS_and_NP([smiles]) and pull out SAS and NP from its table."""
    s = calculate_SAS_and_NP([smiles])
    sas = np_score = None
    for line in s.splitlines():
        if 'SAS Score' in line:
            continue
        # Skip the `|---|---|` separator: after stripping pipes+spaces only '-'/':'
        # should remain.
        if set(re.sub(r'[|\s]', '', line)) <= set('-:'):
            continue
        # Data row: | <smiles> | <sas> | <np> |
        m = re.search(r'([-\d.]+)\s*\|\s*([-\d.]+)\s*\|?\s*$', line)
        if m:
            sas, np_score = float(m.group(1)), float(m.group(2))
    return sas, np_score


# --- sqlite ------------------------------------------------------------------

def open_db(path):
    conn = sqlite3.connect(path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS molecules (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            canonical_smiles TEXT    NOT NULL UNIQUE,
            inchikey         TEXT,
            original_smiles  TEXT,
            source_md        TEXT,
            model            TEXT,
            protein          TEXT,
            docking_score    REAL,
            qed              REAL,
            alogp            REAL,
            np_score         REAL,
            sas_score        REAL,
            added_at         TEXT
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_inchikey ON molecules(inchikey)")
    conn.commit()
    return conn


def is_known(conn, canonical, inchikey):
    """True if either the canonical SMILES or the InChIKey is already in the DB."""
    row = conn.execute(
        "SELECT 1 FROM molecules WHERE canonical_smiles = ? LIMIT 1",
        (canonical,)).fetchone()
    if row:
        return True
    if inchikey:
        row = conn.execute(
            "SELECT 1 FROM molecules WHERE inchikey = ? LIMIT 1",
            (inchikey,)).fetchone()
        if row:
            return True
    return False


# --- main --------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('md_file', nargs='?', help='molopt.py results .md to verify. '
                   'If omitted, the newest results/*.md is used.')
    p.add_argument('--db', default=os.path.join(_HERE, 'molecules.sqlite'),
                   help='sqlite DB path (default: ./molecules.sqlite). Created if absent.')
    p.add_argument('--protein', default=None,
                   help='Docking target. Default: parsed from the md header.')
    p.add_argument('--min-heavy-atoms', type=int, default=5,
                   help='Minimum heavy atoms for a parsed token to count as a proposed '
                        'molecule (filters substituent-fragment noise). Use --no-filter '
                        'to accept any valid SMILES. Default: 5.')
    p.add_argument('--no-filter', action='store_true',
                   help='Accept every valid SMILES (min heavy atoms = 1).')
    p.add_argument('--dry-run', action='store_true',
                   help='List the extracted SMILES without docking or touching the DB.')
    args = p.parse_args(argv)

    md_file = args.md_file
    if not md_file:
        cands = sorted(glob.glob(os.path.join(_HERE, 'results', '*.md')))
        if not cands:
            print("No .md found in results/ and no md_file given.", file=sys.stderr)
            return 2
        md_file = cands[-1]
        print(f"No md_file given; using newest: {md_file}")

    md_file = os.path.abspath(md_file)
    if not os.path.isfile(md_file):
        print(f"Not a file: {md_file}", file=sys.stderr)
        return 2

    with open(md_file, 'r') as f:
        md_text = f.read()

    protein, model = parse_header(md_text)
    protein = args.protein or protein or 'HMGCR'
    if model is None:
        model = os.path.basename(md_file).split('_HMGCR_')[0]  # best-effort fallback

    block = last_model_response_block(md_text)
    if not block:
        print(f"No model-response block found in {md_file} "
              f"(run may have ended before any model turn). Nothing to verify.")
        return 0

    min_hac = 1 if args.no_filter else args.min_heavy_atoms
    molecules = extract_smiles(block, min_hac)

    print(f"Run:       {os.path.basename(md_file)}")
    print(f"Model:      {model}")
    print(f"Protein:    {protein}  (docking target)")
    print(f"Last model-response block: {len(block)} chars")
    print(f"Extracted {len(molecules)} valid SMILES "
          f"(min heavy atoms = {min_hac})")

    if args.dry_run:
        for orig, canon, mol in molecules:
            inchikey = Chem.MolToInchiKey(mol)
            print(f"  {canon}  (HAC={mol.GetNumHeavyAtoms()}, InChIKey={inchikey})  "
                  f"src: {orig!r}")
        return 0

    # Mirror molopt.py: point the shared scoring_args at this run's protein.
    scoring_args[0] = os.cpu_count()
    scoring_args[1] = protein

    conn = open_db(args.db)
    inserted = skipped = failed = 0

    print(f"\nVerifying {len(molecules)} molecules against {protein} ...")
    print("-" * 88)
    print(f"{'canonical SMILES':<50} {'dock':>7} {'QED':>5} {'aLogP':>6} {'SAS':>5} {'NP':>6}  status")
    print("-" * 88)

    for orig, canon, mol in molecules:
        inchikey = Chem.MolToInchiKey(mol)
        if is_known(conn, canon, inchikey):
            skipped += 1
            print(f"{canon[:50]:<50} {'--':>7} {'--':>5} {'--':>6} {'--':>5} {'--':>6}  "
                  f"already in DB, skipped")
            continue

        # Docking (reloads the target each call, as in the live runs).
        try:
            with _quiet_stdout():
                docking_score, _aux = scoring_function(orig)
        except Exception as err:
            docking_score = None
            print(f"[docking failed for {canon}: {err}]", file=sys.stderr)
        # QED + aLogP
        try:
            with _quiet_stdout():
                qed, alogp = _parse_lipinski(orig)
        except Exception:
            qed = alogp = None
        # SAS + NP
        try:
            with _quiet_stdout():
                sas, np_score = _parse_sas_np(orig)
        except Exception:
            sas = np_score = None

        if docking_score is None and qed is None and sas is None:
            failed += 1
            status = "all metrics failed"
        else:
            inserted += 1
            status = "inserted"

        conn.execute("""
            INSERT OR IGNORE INTO molecules
              (canonical_smiles, inchikey, original_smiles, source_md, model, protein,
               docking_score, qed, alogp, np_score, sas_score, added_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """, (canon, inchikey, orig, os.path.basename(md_file), model, protein,
              docking_score, qed, alogp, np_score, sas, datetime.now().isoformat(timespec='seconds')))
        conn.commit()

        d = f"{docking_score:.2f}" if docking_score is not None else "NA"
        q = f"{qed:.3f}" if qed is not None else "NA"
        a = f"{alogp:.2f}" if alogp is not None else "NA"
        s2 = f"{sas:.2f}" if sas is not None else "NA"
        n = f"{np_score:.2f}" if np_score is not None else "NA"
        print(f"{canon[:50]:<50} {d:>7} {q:>5} {a:>6} {s2:>5} {n:>6}  {status}")

    print("-" * 88)
    print(f"\nDone. inserted={inserted}  skipped(already in DB)={skipped}  "
          f"failed={failed}  -> {args.db}")
    conn.close()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())