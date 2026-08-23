#!/usr/bin/env python3
"""
analyze_replicates.py - Compare final compounds across adversary sets / replicates.

Reads a batch manifest written by run_replicates.py, extracts the final-compound SMILES from
each completed run's results .md (reusing verify_results.py's extractor), recomputes the five
metrics (docking, QED, aLogP, SAS, NP) by reusing the project helpers, and writes:

  analysis/compounds_<batch>.csv           per-compound (one row per proposed molecule)
  analysis/summary_<batch>.csv              per-set aggregate stats
  analysis/best_per_replicate_<batch>.csv   best-by-docking molecule per (set, replicate)
  analysis/dock_dist_by_set.png             docking-score distribution per set (box)
  analysis/best_dock_by_replicate.png       min docking score per replicate, by set (strip)
  analysis/qed_vs_dock.png                   QED vs docking scatter, coloured by set
  analysis/property_dist_by_set.png          SAS / NP distributions per set (2-panel violin)

Unlike verify_results.py (which dedups globally and skips known molecules), this calls
extract_smiles per-run independently, so the same SMILES appearing in multiple replicates is
KEPT -- that overlap is part of the comparison.

--skip-docking omits the CPU-heavy dockstring call (leaves the docking column blank) so the
RDKit-only metrics (QED/aLogP/SAS/NP) plus all parsing/CSV/plots can be checked without loading
the CPU. Run the venv activated (source fao-env/bin/activate) for the helpers to import.

Usage:
  python analyze_replicates.py --batch-dir results/batches/<batch_id>
  python analyze_replicates.py --manifest results/batches/<batch_id>/manifest.json --skip-docking
"""

import os
import sys
import json
import glob
import argparse

# --- Path + NumPy-2 shim (mirror molopt.py / verify_results.py) ------------
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, 'code'))

import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

# Reuse the existing extractor + metric parsers (do not reimplement).
from verify_results import (  # noqa: E402
    all_model_response_blocks,
    last_model_response_block,
    extract_smiles,
    _parse_lipinski,
    _parse_sas_np,
    _quiet_stdout,
)
from docking_module import scoring_function, scoring_args  # noqa: E402
from docking_module import contacted_residues, target_residues  # noqa: E402

import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')  # headless: write PNGs, no display needed
import matplotlib.pyplot as plt  # noqa: E402

# The known binding-site residues a proposed molecule should contact to count
# as 'in the pocket', as named in the system message (task_specific_tools).
# Empty if the prompt doesn't state them -- then we fall back to the pose-placed
# signal (aux is not None) for docked_in_pocket.
_TARGET_RESIDUES = target_residues()


# --- set label -> (proposer provider, adversary provider) -------------------

def _providers(set_label):
    """'openai_vs_anthropic' -> ('openai', 'anthropic'); etc."""
    parts = set_label.split('_vs_')
    if len(parts) != 2:
        return (None, None)
    return parts[0], parts[1]


def _model_for(models, provider):
    """Resolve a model name for a provider from the manifest's models dict."""
    if provider is None or not models:
        return None
    return models.get(provider)


# --- per-run extraction -----------------------------------------------------

def extract_run_compounds(md_path, protein, min_heavy_atoms, skip_docking):
    """Parse one run's .md and return a list of per-compound dicts (metrics included).

    Sets the shared scoring_args to this run's protein (mirrors verify_results.py), then
    reuses its block extractor + SMILES extractor + metric parsers.

    The proposer's *final* turn is preferred, but if it produced no parseable SMILES
    (e.g. a proposer cut off mid-thought by the turn cap), fall back to the most recent
    earlier proposer turn that did. Each row carries a `source_turn` (0 = initial
    response, 1..N = adversary-refinement turns) so a fallback is visible in the CSV.
    """
    with open(md_path, 'r') as f:
        md_text = f.read()
    blocks = all_model_response_blocks(md_text)
    if not blocks:
        return []

    # Walk proposer blocks newest -> oldest; use the first that yields >=1 compound.
    chosen, source_turn = None, None
    for idx in range(len(blocks) - 1, -1, -1):
        comps = list(extract_smiles(blocks[idx], min_heavy_atoms))
        if comps:
            chosen, source_turn = comps, idx
            break
    if chosen is None:
        return []

    scoring_args[0] = os.cpu_count()
    scoring_args[1] = protein

    rows = []
    for orig, canon, mol in chosen:
        inchikey = Chem.MolToInchiKey(mol)
        # Docking (CPU-heavy; skippable). scoring_function returns (score, aux);
        # aux is None when dockstring could not place a pose in the pocket box
        # (the same 'Docking failed' check dock_and_get_interacting_residues uses).
        # docked_in_pocket is True when the pose contacts >=1 of the target
        # binding-site residues named in the system message (_TARGET_RESIDUES);
        # if no target residues are defined or contact analysis can't run, it
        # falls back to the pose-placed signal (aux is not None). A failed dock
        # leaves docking as None (not a misleading 0.0).
        docking = None
        docked_in_pocket = None
        n_target_contacts = None
        if not skip_docking:
            try:
                with _quiet_stdout():
                    score, aux = scoring_function(orig)
                docking = score if aux is not None else None
                if aux is None:
                    docked_in_pocket = False
                elif _TARGET_RESIDUES:
                    contacts = contacted_residues(aux)  # set, or None if unavailable
                    if contacts is None:
                        # No receptor PDB / contact analysis unavailable: pose was
                        # placed, but we can't check which residues it contacted.
                        docked_in_pocket = True
                    else:
                        n_target_contacts = len(contacts & _TARGET_RESIDUES)
                        docked_in_pocket = n_target_contacts > 0
                else:
                    # No target residues defined for this protein -> pose-placed signal.
                    docked_in_pocket = True
            except Exception:
                docking = None
                docked_in_pocket = False
                n_target_contacts = None
        # QED + aLogP (RDKit-only, cheap).
        try:
            with _quiet_stdout():
                qed, alogp = _parse_lipinski(orig)
        except Exception:
            qed = alogp = None
        # SAS + NP (RDKit-only, cheap).
        try:
            with _quiet_stdout():
                sas, np_score = _parse_sas_np(orig)
        except Exception:
            sas = np_score = None
        rows.append({
            'original_smiles': orig,
            'canonical_smiles': canon,
            'inchikey': inchikey,
            'docking': docking,
            'docked_in_pocket': docked_in_pocket,
            'n_target_contacts': n_target_contacts,
            'qed': qed,
            'alogp': alogp,
            'sas': sas,
            'np': np_score,
            'source_turn': source_turn,
        })
    return rows


# --- per-run tool-call usage (from the JSON sidecar) ------------------------

# The tool that actually runs Vina (the CPU-costly one); count it separately.
_DOCK_TOOL = 'dock_and_get_interacting_residues'


def _msg_tool_calls(msg):
    """Tool names called in one message, in either native format.

    Anthropic: content is a list of blocks, tool calls are {type:'tool_use',...}.
    OpenAI:    tool calls are a top-level 'tool_calls' list. Either may be absent.
    Gemini:   'parts' is a list; a tool call is a part dict with a 'function_call'
              dict (carrying 'name'), optionally with a 'thought_signature'.
    """
    names = []
    # OpenAI native format.
    for tc in (msg.get('tool_calls') or []):
        try:
            names.append(tc['function']['name'])
        except Exception:
            pass
    # Anthropic native format (content blocks).
    content = msg.get('content')
    if isinstance(content, list):
        for b in content:
            if isinstance(b, dict) and b.get('type') == 'tool_use':
                names.append(b.get('name') or '')
    # Gemini native format (parts list with function_call dicts).
    parts = msg.get('parts')
    if isinstance(parts, list):
        for b in parts:
            if isinstance(b, dict) and isinstance(b.get('function_call'), dict):
                names.append(b['function_call'].get('name') or '')
    return names


def _is_tool_result_user(msg):
    """Is this 'user' message a tool-result/over-cap stub (intra-turn), not a
    real prompt that starts a new proposer turn? Handles all native formats."""
    c = msg.get('content')
    if isinstance(c, list):
        return any(isinstance(b, dict) and b.get('type') == 'tool_result' for b in c)
    if isinstance(c, str):
        # OpenAI tool results are separate 'tool' role msgs (not 'user'), but an
        # over-cap nudge is a 'user' string starting with the cap message.
        return c.startswith('Tool-call limit') or c.startswith('You have reached the tool')
    # Gemini native format: a 'user' message whose parts carry a function_response
    # is an intra-turn tool result, not a real proposer-turn prompt.
    parts = msg.get('parts')
    if isinstance(parts, list):
        return any(isinstance(b, dict) and b.get('function_response') is not None for b in parts)
    return False


def tool_usage_from_sidecar(sidecar_path):
    """Summarise tool usage for one run from its JSON sidecar.

    Returns a dict {n_turns, total_tool_rounds, total_tool_calls, n_docks,
    rounds_per_turn, calls_per_turn} or None if the sidecar can't be read.
    A 'tool round' = one assistant message that carried >=1 tool call (each
    such message is a separate API round-trip). n_turns = number of proposer
    turns = count of non-tool-result user messages (the seed prompt + each
    adversary-feedback prompt), which is written_at_turn + 1 (the sidecar's
    written_at_turn excludes the initial seed turn).
    """
    try:
        with open(sidecar_path) as f:
            d = json.load(f)
    except Exception:
        return None
    msgs = d.get('messages') or []
    n_turns = sum(1 for m in msgs
                  if m.get('role') == 'user' and not _is_tool_result_user(m))
    if not n_turns:  # fall back to the sidecar's own counter
        n_turns = d.get('written_at_turn') or 0
    rounds = calls = docks = 0
    for m in msgs:
        # 'assistant' for OpenAI/Anthropic; 'model' for Gemini's native format.
        if m.get('role') not in ('assistant', 'model'):
            continue
        names = _msg_tool_calls(m)
        if names:
            rounds += 1
            calls += len(names)
            docks += sum(1 for n in names if n == _DOCK_TOOL)
    return {
        'n_turns': n_turns,
        'total_tool_rounds': rounds,
        'total_tool_calls': calls,
        'n_docks': docks,
        'rounds_per_turn': round(rounds / n_turns, 2) if n_turns else None,
        'calls_per_turn': round(calls / n_turns, 2) if n_turns else None,
    }


# --- stats ------------------------------------------------------------------

_METRIC_COLS = ['docking', 'qed', 'alogp', 'sas', 'np']


def summary_stats(df):
    """One row per set_label with aggregate stats across replicates."""
    out = []
    for set_label, g in df.groupby('set_label', sort=True):
        rec = {
            'set_label': set_label,
            'n_replicates': g['replicate'].nunique(),
            'n_compounds': len(g),
            'n_unique': g['canonical_smiles'].nunique(),
        }
        for col in _METRIC_COLS:
            s = g[col].dropna()
            rec[f'{col}_mean'] = round(s.mean(), 3) if len(s) else None
            if col == 'docking':
                rec['docking_median'] = round(s.median(), 3) if len(s) else None
                rec['docking_best'] = round(s.min(), 3) if len(s) else None  # lower = better
        out.append(rec)
    return pd.DataFrame(out).set_index('set_label')


def best_per_replicate(df):
    """Best-by-docking molecule per (set, replicate). Requires docking scores."""
    out = []
    for (set_label, rep), g in df.groupby(['set_label', 'replicate'], sort=True):
        gd = g.dropna(subset=['docking'])
        if len(gd):
            best = gd.loc[gd['docking'].idxmin()]
        else:
            best = g.iloc[0]  # no docking: keep a representative row (docking NaN)
        out.append({
            'set_label': set_label,
            'replicate': rep,
            'canonical_smiles': best['canonical_smiles'],
            'docking': best['docking'],
            'docked_in_pocket': best.get('docked_in_pocket'),
            'n_target_contacts': best.get('n_target_contacts'),
            'qed': best['qed'],
            'alogp': best['alogp'],
            'sas': best['sas'],
            'np': best['np'],
        })
    return pd.DataFrame(out)


# --- plots ------------------------------------------------------------------

def _plot_dock_dist(df, out_path):
    sets = sorted(df['set_label'].dropna().unique())
    data = [df.loc[df['set_label'] == s, 'docking'].dropna().values for s in sets]
    if not any(len(d) for d in data):
        return False
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(data, tick_labels=sets, showfliers=True)
    ax.set_ylabel('docking score (lower = better)')
    ax.set_title('Final-compound docking score by adversary set')
    ax.tick_params(axis='x', rotation=20)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


def _plot_best_by_replicate(df, out_path):
    sub = df.dropna(subset=['docking'])
    if sub.empty:
        return False
    fig, ax = plt.subplots(figsize=(8, 5))
    sets = sorted(sub['set_label'].unique())
    for s in sets:
        g = sub[sub['set_label'] == s].sort_values('replicate')
        ax.plot(g['replicate'], g['docking'], marker='o', linestyle='-', label=s, alpha=0.8)
    ax.set_xlabel('replicate')
    ax.set_ylabel('best docking score (min, lower = better)')
    ax.set_title('Best docking score per replicate, by adversary set')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


def _plot_qed_vs_dock(df, out_path):
    sub = df.dropna(subset=['docking', 'qed'])
    if sub.empty:
        return False
    fig, ax = plt.subplots(figsize=(7, 5))
    sets = sorted(sub['set_label'].unique())
    for s in sets:
        g = sub[sub['set_label'] == s]
        ax.scatter(g['docking'], g['qed'], label=s, alpha=0.6, s=30)
    ax.set_xlabel('docking score (lower = better)')
    ax.set_ylabel('QED (higher = more drug-like)')
    ax.set_title('QED vs docking by adversary set')
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
    fig.suptitle('Final-compound SAS / NP by adversary set')
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return True


# --- main -------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog='analyze_replicates.py',
        description='Compare final compounds across adversary sets / replicates.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--batch-dir', default=None,
                   help='A results/batches/<batch_id> dir (reads its manifest.json).')
    p.add_argument('--manifest', default=None,
                   help='Path to a manifest.json (alternative to --batch-dir).')
    p.add_argument('--out-dir', default=None,
                   help='Where to write CSVs + PNGs (default: <batch-dir>/analysis).')
    p.add_argument('--min-heavy-atoms', type=int, default=5,
                   help='Min heavy atoms for a parsed SMILES to count (passed to extract_smiles). '
                        'Default: 5.')
    p.add_argument('--skip-docking', action='store_true',
                   help='Skip the CPU-heavy dockstring call (leave docking blank). QED/aLogP/'
                        'SAS/NP (RDKit-only) and all CSVs/plots are still produced.')
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
    out_dir = args.out_dir or os.path.join(batch_dir, 'analysis')
    os.makedirs(out_dir, exist_ok=True)

    entries = [e for e in manifest.get('entries', []) if e.get('status') == args.status]
    if not entries:
        print(f"No '{args.status}' entries in manifest. "
              f"Statuses present: {sorted({e.get('status') for e in manifest.get('entries', [])})}")
        return 0

    print(f"Analyzing {len(entries)} run(s) from batch {batch_id} "
          f"{'[skip-docking]' if args.skip_docking else '[with docking]'}")

    rows = []
    tool_rows = []  # one per rep: tool-call usage from the JSON sidecar
    for e in entries:
        md = e.get('md_path')
        if not md or not os.path.isfile(md):
            print(f"  [skip {e['set_label']} rep{e['replicate']}: md_path missing/absent]")
            continue
        protein = e.get('protein', 'HMGCR')
        comps = extract_run_compounds(md, protein, args.min_heavy_atoms, args.skip_docking)
        pp, ap = _providers(e['set_label'])
        for c in comps:
            rows.append({
                'set_label': e['set_label'],
                'replicate': e['replicate'],
                'proposer_provider': pp,
                'proposer_model': _model_for(models, pp),
                'adversary_provider': ap,
                'adversary_model': _model_for(models, ap),
                'protein': protein,
                **c,
            })
        n_turns = len(all_model_response_blocks(open(md).read())) if md else 0
        t = comps[0]['source_turn'] if comps else None
        note = '' if (t is None or n_turns == 0 or t == n_turns - 1) else f' [fallback: proposer turn {t}, last turn had no SMILES]'
        n_pocket = sum(1 for c in comps if c.get('docked_in_pocket'))
        pocket_note = '' if args.skip_docking else f' ({n_pocket}/{len(comps)} in pocket)'
        # Tool-call usage from the sibling sidecar (independent of docking).
        tu = tool_usage_from_sidecar(e.get('sidecar_path'))
        tool_note = ''
        if tu is not None:
            tool_rows.append({
                'set_label': e['set_label'],
                'replicate': e['replicate'],
                'proposer_model': _model_for(models, pp),
                **tu,
            })
            tool_note = (f" | tools: {tu['total_tool_rounds']} rounds / "
                         f"{tu['total_tool_calls']} calls ({tu['n_docks']} docks) "
                         f"over {tu['n_turns']} turns "
                         f"= {tu['rounds_per_turn']}/turn")
        print(f"  {e['set_label']} rep{e['replicate']}: {len(comps)} compounds{pocket_note}{note}{tool_note}")

    if not rows:
        print("No compounds extracted from any run. Nothing to write.")
        return 0

    df = pd.DataFrame(rows)
    stem = batch_id
    compounds_csv = os.path.join(out_dir, f'compounds_{stem}.csv')
    summary_csv = os.path.join(out_dir, f'summary_{stem}.csv')
    best_csv = os.path.join(out_dir, f'best_per_replicate_{stem}.csv')

    df.to_csv(compounds_csv, index=False)
    print(f"\nWrote per-compound CSV: {compounds_csv} ({len(df)} rows)")

    summ = summary_stats(df)
    summ.to_csv(summary_csv)
    print(f"Wrote per-set summary:  {summary_csv}")
    print(summ.to_string())

    if not args.skip_docking:
        best = best_per_replicate(df)
        best.to_csv(best_csv, index=False)
        print(f"Wrote best-per-replicate: {best_csv}")

    # Per-replicate tool-call usage (always emitted; sourced from sidecars).
    if tool_rows:
        tool_csv = os.path.join(out_dir, f'tool_usage_{stem}.csv')
        pd.DataFrame(tool_rows).to_csv(tool_csv, index=False)
        print(f"Wrote tool-usage per replicate: {tool_csv}")

    # Plots.
    made = []
    if _plot_dock_dist(df, os.path.join(out_dir, 'dock_dist_by_set.png')):
        made.append('dock_dist_by_set.png')
    if not args.skip_docking and _plot_best_by_replicate(df, os.path.join(out_dir, 'best_dock_by_replicate.png')):
        made.append('best_dock_by_replicate.png')
    if _plot_qed_vs_dock(df, os.path.join(out_dir, 'qed_vs_dock.png')):
        made.append('qed_vs_dock.png')
    if _plot_property_dist(df, os.path.join(out_dir, 'property_dist_by_set.png')):
        made.append('property_dist_by_set.png')
    if made:
        print(f"Wrote plots: {', '.join(made)}")
    else:
        print("No plots produced (insufficient data).")

    print(f"\nDone. Outputs in: {out_dir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())