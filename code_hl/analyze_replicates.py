#!/usr/bin/env python3
"""
analyze_replicates.py - Compare final compounds across sets / replicates (HL-gap version).

The HOMO-LUMO counterpart of code/analyze_replicates.py. Reads a batch manifest written by
code_hl/run_replicates.py (agentic sets) or code_hl/run_zero_few_shot.py (zero/few/frag-shot
baselines), extracts each run's final-compound SMILES, recomputes the metrics (HL gap, QED,
aLogP, SAS, NP) and writes:

  analysis/compounds_<batch>.csv            per-compound (one row per proposed molecule)
  analysis/summary_<batch>.csv              per-set aggregate stats
  analysis/best_per_replicate_<batch>.csv   best-by-gap molecule per (set, replicate)
  analysis/tool_usage_<batch>.csv           per-replicate tool-call usage (agentic runs)
  analysis/gap_dist_by_set.png              gap distribution per set (box)
  analysis/best_gap_by_replicate.png        min gap per replicate, by set
  analysis/qed_vs_gap.png                   QED vs gap scatter, coloured by set
  analysis/property_dist_by_set.png         SAS / NP distributions per set (2-panel violin)

Differences from the docking version, all forced by the objective:
  - the metric is the GFN2-xTB HOMO-LUMO gap in eV, recomputed with
    hl_gap_module.scoring_function. Lower is better, same direction as docking.
  - scoring_function returns its 100.0 failure sentinel (invalid SMILES, 3D embed
    failure, non-convergence) rather than raising; those become gap=None, so a failed
    calculation is never averaged in as a real (terrible) value.
  - there is no pocket: the docking version's docked_in_pocket / n_target_contacts
    columns have no analogue and are dropped.
  - the expensive tool counted per replicate is calculate_HL_gap, not the docking tool.

The SMILES extractor and property parsers are imported from code/verify_results.py rather
than reimplemented, so HL and docking runs are parsed by identical code and stay comparable.
That module imports docking_module at module level, so dockstring/oddt load here too even
though nothing docks -- harmless, and cheaper than letting the two arms' parsers drift.

Gap recomputation is ~0.02-0.1 s/molecule, so unlike the docking version there is no need to
skip it; --skip-gaps exists only for parser/CSV/plot checks.

Usage:
  python analyze_replicates.py --batch-dir results/batches/hl_batches/<batch_id>
  python analyze_replicates.py --manifest .../manifest.json --skip-gaps
"""

import os
import sys
import json
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))   # code_hl/
_ROOT = os.path.dirname(_HERE)                        # repo root
# code_hl first: 'MolPropOp' must resolve to the HL copy (its lipinski is identical to the
# docking copy's, but this keeps the HL side self-consistent). 'docking_module', which
# verify_results imports, exists only in code/, so it resolves there.
sys.path.insert(0, os.path.join(_ROOT, 'code'))
sys.path.insert(0, _HERE)

import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin

from rdkit import Chem, RDLogger  # noqa: E402
RDLogger.DisableLog('rdApp.*')

# Reuse the docking arm's extractor + metric parsers (do not reimplement).
from verify_results import (  # noqa: E402
    all_model_response_blocks,
    extract_smiles,
    _parse_lipinski,
    _parse_sas_np,
    _quiet_stdout,
)
from hl_gap_module import scoring_function, scoring_args  # noqa: E402

import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

# scoring_function's failure sentinel; any value at or above this is 'no result',
# never a real gap (real gaps are single-digit eV).
_FAIL_SENTINEL = 100.0

# The tool that actually runs GFN2-xTB (the costly one); counted separately.
_GAP_TOOL = 'calculate_HL_gap'


# --- set label -> (proposer provider, adversary provider) -------------------

def _providers(set_label):
    """'openai_vs_anthropic' -> ('openai', 'anthropic').

    Baseline manifests use a bare model label ('anthropic', 'kimi-k2.6') with no
    adversary, so there is nothing to split: the proposer is the label itself.
    """
    parts = set_label.split('_vs_')
    if len(parts) != 2:
        return (set_label, None)
    return parts[0], parts[1]


def _model_for(models, provider):
    if provider is None or not models:
        return None
    return models.get(provider)


# --- per-run extraction -----------------------------------------------------

def extract_run_compounds(md_path, method, min_heavy_atoms, skip_gaps):
    """Parse one run's .md and return a list of per-compound dicts with metrics.

    The proposer's *final* turn is preferred; if it yielded no parseable SMILES (a
    proposer cut off by the turn cap, or one that answered with questions instead of
    molecules), fall back to the most recent earlier turn that did. `source_turn`
    records which, so a fallback is visible in the CSV.
    """
    with open(md_path, 'r') as f:
        md_text = f.read()
    blocks = all_model_response_blocks(md_text)
    if not blocks:
        return []

    chosen, source_turn = None, None
    for idx in range(len(blocks) - 1, -1, -1):
        comps = list(extract_smiles(blocks[idx], min_heavy_atoms))
        if comps:
            chosen, source_turn = comps, idx
            break
    if chosen is None:
        return []

    scoring_args[0] = os.cpu_count()
    scoring_args[1] = method

    rows = []
    for orig, canon, mol in chosen:
        inchikey = Chem.MolToInchiKey(mol)
        gap = None
        if not skip_gaps:
            try:
                with _quiet_stdout():
                    score, _aux = scoring_function(orig)
                gap = None if (score is None or score >= _FAIL_SENTINEL) else float(score)
            except Exception:
                gap = None
        try:
            with _quiet_stdout():
                qed, alogp = _parse_lipinski(orig)
        except Exception:
            qed = alogp = None
        try:
            with _quiet_stdout():
                sas, np_score = _parse_sas_np(orig)
        except Exception:
            sas = np_score = None
        rows.append({
            'original_smiles': orig,
            'canonical_smiles': canon,
            'inchikey': inchikey,
            'gap': gap,
            'qed': qed,
            'alogp': alogp,
            'sas': sas,
            'np': np_score,
            'source_turn': source_turn,
        })
    return rows


# --- per-run tool-call usage (from the JSON sidecar) ------------------------

def _msg_tool_calls(msg):
    """Tool names called in one message, in any of the native formats.

    OpenAI/Ollama: top-level 'tool_calls' list. Anthropic: content blocks of
    type 'tool_use'. Gemini: 'parts' entries carrying a 'function_call' dict.
    """
    names = []
    for tc in (msg.get('tool_calls') or []):
        try:
            names.append(tc['function']['name'])
        except Exception:
            pass
    content = msg.get('content')
    if isinstance(content, list):
        for b in content:
            if isinstance(b, dict) and b.get('type') == 'tool_use':
                names.append(b.get('name') or '')
    parts = msg.get('parts')
    if isinstance(parts, list):
        for b in parts:
            if isinstance(b, dict) and isinstance(b.get('function_call'), dict):
                names.append(b['function_call'].get('name') or '')
    return names


def _is_tool_result_user(msg):
    """Is this 'user' message an intra-turn tool-result/nudge stub rather than a
    real prompt starting a new proposer turn?"""
    c = msg.get('content')
    if isinstance(c, list):
        return any(isinstance(b, dict) and b.get('type') == 'tool_result' for b in c)
    if isinstance(c, str):
        return (c.startswith('Tool-call limit') or c.startswith('You have reached the tool')
                or c.startswith('The previous call failed with a connection error'))
    parts = msg.get('parts')
    if isinstance(parts, list):
        return any(isinstance(b, dict) and b.get('function_response') is not None for b in parts)
    return False


def tool_usage_from_sidecar(sidecar_path):
    """Summarise one run's tool usage. n_turns counts refinement turns only, so it is
    directly comparable to --max-turns. Returns None if the sidecar can't be read."""
    try:
        with open(sidecar_path) as f:
            d = json.load(f)
    except Exception:
        return None
    msgs = d.get('messages') or []
    n_proposer = sum(1 for m in msgs
                     if m.get('role') == 'user' and not _is_tool_result_user(m))
    n_turns = max(n_proposer - 1, 0)
    if not n_proposer:
        n_turns = d.get('written_at_turn') or 0
    rounds = calls = gap_calls = 0
    for m in msgs:
        if m.get('role') not in ('assistant', 'model'):
            continue
        names = _msg_tool_calls(m)
        if names:
            rounds += 1
            calls += len(names)
            gap_calls += sum(1 for n in names if n == _GAP_TOOL)
    return {
        'n_turns': n_turns,
        'total_tool_rounds': rounds,
        'total_tool_calls': calls,
        'n_gap_calls': gap_calls,
        'rounds_per_turn': round(rounds / n_proposer, 2) if n_proposer else None,
        'calls_per_turn': round(calls / n_proposer, 2) if n_proposer else None,
    }


# --- stats ------------------------------------------------------------------

_METRIC_COLS = ['gap', 'qed', 'alogp', 'sas', 'np']


def summary_stats(df, skip_gaps=False):
    """One row per set_label with aggregate stats across replicates.

    n_gap_failed is None under --skip-gaps: nothing was computed, so an all-NaN gap
    column means 'not measured', not 'every calculation failed'.
    """
    out = []
    for set_label, g in df.groupby('set_label', sort=True):
        rec = {
            'set_label': set_label,
            'n_replicates': g['replicate'].nunique(),
            'n_compounds': len(g),
            'n_unique': g['canonical_smiles'].nunique(),
            'n_gap_failed': None if skip_gaps else int(g['gap'].isna().sum()),
        }
        for col in _METRIC_COLS:
            s = g[col].dropna()
            rec[f'{col}_mean'] = round(s.mean(), 3) if len(s) else None
            if col == 'gap':
                rec['gap_median'] = round(s.median(), 3) if len(s) else None
                rec['gap_best'] = round(s.min(), 3) if len(s) else None  # lower = better
        out.append(rec)
    return pd.DataFrame(out).set_index('set_label')


def best_per_replicate(df):
    """Smallest-gap molecule per (set, replicate)."""
    out = []
    for (set_label, rep), g in df.groupby(['set_label', 'replicate'], sort=True):
        gd = g.dropna(subset=['gap'])
        best = gd.loc[gd['gap'].idxmin()] if len(gd) else g.iloc[0]
        out.append({
            'set_label': set_label,
            'replicate': rep,
            'canonical_smiles': best['canonical_smiles'],
            'gap': best['gap'],
            'qed': best['qed'],
            'alogp': best['alogp'],
            'sas': best['sas'],
            'np': best['np'],
        })
    return pd.DataFrame(out)


# --- plots ------------------------------------------------------------------

def _plot_gap_dist(df, out_path):
    sets = sorted(df['set_label'].dropna().unique())
    data = [df.loc[df['set_label'] == s, 'gap'].dropna().values for s in sets]
    if not any(len(d) for d in data):
        return False
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(data, tick_labels=sets, showfliers=True)
    ax.set_ylabel('HOMO-LUMO gap / eV (lower = better)')
    ax.set_title('Final-compound HOMO-LUMO gap by set')
    ax.tick_params(axis='x', rotation=20)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


def _plot_best_by_replicate(df, out_path):
    sub = df.dropna(subset=['gap'])
    if sub.empty:
        return False
    fig, ax = plt.subplots(figsize=(8, 5))
    for s in sorted(sub['set_label'].unique()):
        g = sub[sub['set_label'] == s].groupby('replicate')['gap'].min().reset_index()
        ax.plot(g['replicate'], g['gap'], marker='o', linestyle='-', label=s, alpha=0.8)
    ax.set_xlabel('replicate')
    ax.set_ylabel('best gap / eV (min, lower = better)')
    ax.set_title('Best HOMO-LUMO gap per replicate, by set')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


def _plot_qed_vs_gap(df, out_path):
    sub = df.dropna(subset=['gap', 'qed'])
    if sub.empty:
        return False
    fig, ax = plt.subplots(figsize=(7, 5))
    for s in sorted(sub['set_label'].unique()):
        g = sub[sub['set_label'] == s]
        ax.scatter(g['gap'], g['qed'], label=s, alpha=0.6, s=30)
    ax.set_xlabel('HOMO-LUMO gap / eV (lower = better)')
    ax.set_ylabel('QED')
    ax.set_title('QED vs HOMO-LUMO gap by set')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


def _plot_property_dist(df, out_path):
    sub = df.dropna(subset=['sas', 'np'])
    if sub.empty:
        return False
    sets = sorted(df['set_label'].unique())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
    for ax, col, ylabel in [(a1, 'sas', 'SAS (1 easy - 10 hard)'),
                            (a2, 'np', 'NP (-5..5, higher = more NP-like)')]:
        data = [df.loc[df['set_label'] == s, col].dropna().values for s in sets]
        if any(len(d) for d in data):
            ax.violinplot(data, showmeans=False, showmedians=True)
            ax.set_xticks(range(1, len(sets) + 1))
            ax.set_xticklabels(sets, rotation=20)
            ax.set_ylabel(ylabel)
    fig.suptitle('Final-compound SAS / NP by set')
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


# --- main -------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog='analyze_replicates.py',
        description='Compare final compounds across sets / replicates (HOMO-LUMO gap).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--batch-dir', default=None,
                   help='A results/batches/hl_batches/<batch_id> dir (reads its manifest.json).')
    p.add_argument('--manifest', default=None,
                   help='Path to a manifest.json (alternative to --batch-dir).')
    p.add_argument('--out-dir', default=None,
                   help='Where to write CSVs + PNGs (default: <batch-dir>/analysis).')
    p.add_argument('--min-heavy-atoms', type=int, default=5,
                   help='Min heavy atoms for a parsed SMILES to count (default: 5).')
    p.add_argument('--skip-gaps', action='store_true',
                   help='Skip the GFN2-xTB recomputation (leave gap blank); parser/CSV/plot check.')
    p.add_argument('--status', default='complete',
                   help='Manifest entry status to include (default: complete).')
    args = p.parse_args(argv)

    manifest_path = args.manifest
    if manifest_path is None:
        if not args.batch_dir:
            print("Give --batch-dir or --manifest.", file=sys.stderr)
            return 2
        manifest_path = os.path.join(args.batch_dir, 'manifest.json')
    if not os.path.isfile(manifest_path):
        print(f"Manifest not found: {manifest_path}", file=sys.stderr)
        return 2
    with open(manifest_path) as f:
        manifest = json.load(f)
    batch_id = manifest.get('batch_id', 'batch')
    batch_dir = manifest.get('batch_dir', os.path.dirname(manifest_path))
    models = manifest.get('models', {})
    default_method = manifest.get('method', 'GFN2-xTB')
    out_dir = args.out_dir or os.path.join(batch_dir, 'analysis')
    os.makedirs(out_dir, exist_ok=True)

    # Skipped entries are complete runs the resume logic did not re-run, so they carry
    # real results and belong in the analysis alongside 'complete'.
    wanted = {args.status, 'skipped'} if args.status == 'complete' else {args.status}
    entries = [e for e in manifest.get('entries', []) if e.get('status') in wanted]
    if not entries:
        print(f"No {sorted(wanted)} entries in manifest. "
              f"Statuses present: {sorted({e.get('status') for e in manifest.get('entries', [])})}")
        return 0

    print(f"Analyzing {len(entries)} run(s) from batch {batch_id} "
          f"{'[skip-gaps]' if args.skip_gaps else '[with GFN2-xTB]'}")

    rows, tool_rows = [], []
    stem = batch_id
    compounds_csv = os.path.join(out_dir, f'compounds_{stem}.csv')
    for e in entries:
        md = e.get('md_path')
        if not md or not os.path.isfile(md):
            print(f"  [skip {e.get('set_label')} rep{e.get('replicate')}: md_path missing]")
            continue
        method = e.get('method', default_method)
        comps = extract_run_compounds(md, method, args.min_heavy_atoms, args.skip_gaps)
        pp, ap = _providers(e['set_label'])
        for c in comps:
            rows.append({
                'set_label': e['set_label'],
                'replicate': e['replicate'],
                'proposer_provider': pp,
                'proposer_model': _model_for(models, pp),
                'adversary_provider': ap,
                'adversary_model': _model_for(models, ap),
                'method': method,
                **c,
            })
        n_blocks = len(all_model_response_blocks(open(md).read()))
        t = comps[0]['source_turn'] if comps else None
        note = '' if (t is None or n_blocks == 0 or t == n_blocks - 1) else \
            f' [fallback: turn {t}, last turn had no SMILES]'
        n_failed = sum(1 for c in comps if c['gap'] is None)
        gap_note = '' if args.skip_gaps else f' ({n_failed} gap failures)' if n_failed else ''
        tu = tool_usage_from_sidecar(e.get('sidecar_path'))
        tool_note = ''
        if tu is not None:
            tool_rows.append({
                'set_label': e['set_label'], 'replicate': e['replicate'],
                'proposer_model': _model_for(models, pp), **tu,
            })
            tool_note = (f" | tools: {tu['total_tool_rounds']} rounds / {tu['total_tool_calls']} calls "
                         f"({tu['n_gap_calls']} gap) over {tu['n_turns']} turns")
        print(f"  {e['set_label']} rep{e['replicate']}: {len(comps)} compounds{gap_note}{note}{tool_note}")

        # Checkpoint after each entry so a killed process keeps what it computed.
        if rows:
            pd.DataFrame(rows).to_csv(compounds_csv, index=False)

    if not rows:
        print("No compounds extracted from any run. Nothing to write.")
        return 0

    df = pd.DataFrame(rows)
    df.to_csv(compounds_csv, index=False)
    print(f"\nWrote per-compound CSV: {compounds_csv} ({len(df)} rows)")

    summ = summary_stats(df, skip_gaps=args.skip_gaps)
    summary_csv = os.path.join(out_dir, f'summary_{stem}.csv')
    summ.to_csv(summary_csv)
    print(f"Wrote per-set summary:  {summary_csv}")
    print(summ.to_string())

    if not args.skip_gaps:
        best_csv = os.path.join(out_dir, f'best_per_replicate_{stem}.csv')
        best_per_replicate(df).to_csv(best_csv, index=False)
        print(f"Wrote best-per-replicate: {best_csv}")

    if tool_rows:
        tool_csv = os.path.join(out_dir, f'tool_usage_{stem}.csv')
        pd.DataFrame(tool_rows).to_csv(tool_csv, index=False)
        print(f"Wrote tool-usage per replicate: {tool_csv}")

    made = []
    if _plot_gap_dist(df, os.path.join(out_dir, 'gap_dist_by_set.png')):
        made.append('gap_dist_by_set.png')
    if not args.skip_gaps and _plot_best_by_replicate(df, os.path.join(out_dir, 'best_gap_by_replicate.png')):
        made.append('best_gap_by_replicate.png')
    if _plot_qed_vs_gap(df, os.path.join(out_dir, 'qed_vs_gap.png')):
        made.append('qed_vs_gap.png')
    if _plot_property_dist(df, os.path.join(out_dir, 'property_dist_by_set.png')):
        made.append('property_dist_by_set.png')
    print(f"Wrote plots: {', '.join(made)}" if made else "No plots produced (insufficient data).")

    print(f"\nDone. Outputs in: {out_dir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
