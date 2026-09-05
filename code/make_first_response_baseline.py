#!/usr/bin/env python3
"""make_first_response_baseline.py - "First-response" baseline for both study arms.

A new baseline condition: what each agentic proposer produces after exactly ONE
response — its tools in the first turn, but no adversary critique and no
multi-turn refinement. Instead of re-running the loop, this mines the existing
agentic transcripts: for each replicate it takes the proposer's FIRST response
— every model text segment in turn 1, i.e. everything in the transcript before
the first "# Adversary feedback" section (turn 1 may interleave text with
intra-turn tool calls) — extracts the proposed SMILES with the study's
canonical extractor (verify_results.extract_smiles, min 5 heavy atoms, deduped
by canonical SMILES), and scores them with that arm's real engine — exactly the
way analyze_replicates.py scores final-turn compounds:

  --study dock  docking_module.scoring_function (Vina via dockstring), plus the
                docked_in_pocket / n_target_contacts pocket check. CPU-heavy
                (~1 dock per compound): run detached, the script is resumable.
  --study hl    hl_gap_module.scoring_function (GFN2-xTB via tblite); failed
                calculations leave the gap blank, per study convention.

Sources are the five primary agentic 5x4 batches per arm. The first reply never
sees an adversary, so for docking the cross-critic batches are used (OpenAI's
row is critiqued by Anthropic there — irrelevant for turn 1); for HL the
self-critic batches. Output is a single baseline-style CSV per arm (one row per
compound, set_label = proposer key, source_turn = 0):

  results/batches/first_response_5x4/analysis/compounds_first_response_5x4.csv
  results/batches/hl_batches/first_response_5x4/analysis/compounds_first_response_5x4.csv
plus summary_<...>.csv and best_per_replicate_<...>.csv in the same schemas as
the per-batch analysis CSVs. Resumable: (set_label, replicate, canonical_smiles)
rows already in the CSV are skipped, so an interrupted docking run loses nothing.

Usage:
  fao-env/bin/python code/make_first_response_baseline.py --study hl
  fao-env/bin/python code/make_first_response_baseline.py --study dock
"""
import argparse
import contextlib
import csv
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)
sys.path.append(_HERE)
sys.path.append(os.path.join(_HERE, '..', 'code_hl'))

from verify_results import all_model_response_blocks, extract_smiles  # noqa: E402

_BATCHES = os.path.join(_ROOT, 'results', 'batches')

# batch_id -> (proposer set_label / provider key, proposer_model, study)
DOCK_BATCHES = {
    'openai_gpt-5.2_vs_anthropic_5x4': ('openai', 'gpt-5.2'),
    'anthropic-haiku-4-5_vs_openai_5x4': ('anthropic', 'claude-haiku-4-5-20251001'),
    'gemini-3-flash-preview_vs_openai_5x4': ('gemini', 'gemini-3-flash-preview'),
    'ollama_kimi-k2.6_vs_openai_5x4': ('kimi', 'kimi-k2.6'),
    'ollama_deepseek-v4-pro_vs_openai_5x4': ('deepseek', 'deepseek-v4-pro'),
}
HL_BATCHES = {
    'hl_batches/hl_gpt-5.2_vs_gpt-5.2_5x4': ('openai', 'gpt-5.2'),
    'hl_batches/hl_claude-haiku-4-5_vs_claude-haiku-4-5_5x4': ('anthropic', 'claude-haiku-4-5-20251001'),
    'hl_batches/hl_gemini-3-flash-preview_vs_gemini_5x4': ('gemini', 'gemini-3-flash-preview'),
    'hl_batches/hl_kimi-k2.6_vs_kimi-k2.6_5x4': ('kimi', 'kimi-k2.6'),
    'hl_batches/hl_deepseek-v4-pro_vs_deepseek-v4-pro_5x4': ('deepseek', 'deepseek-v4-pro'),
}
BATCHES = {'dock': DOCK_BATCHES, 'hl': HL_BATCHES}
OUT_DIR = {'dock': os.path.join(_BATCHES, 'first_response_5x4', 'analysis'),
           'hl': os.path.join(_BATCHES, 'hl_batches', 'first_response_5x4', 'analysis')}
STEM = 'first_response_5x4'
MIN_HEAVY_ATOMS = 5  # same as analyze_replicates.py's default


@contextlib.contextmanager
def _quiet_stdout():
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            yield


def _score_dock(orig):
    """Mirror code/analyze_replicates.py's scoring block for one compound."""
    from docking_module import scoring_function, contacted_residues, target_residues, scoring_args
    scoring_args[0] = os.cpu_count()
    scoring_args[1] = 'HMGCR'
    with _quiet_stdout():
        score, aux = scoring_function(orig)
    docking = score if aux is not None else None
    if aux is None:
        return docking, False, None
    _targets = target_residues()
    if not _targets:
        return docking, True, None
    contacts = contacted_residues(aux)
    if contacts is None:
        return docking, True, None
    n = len(contacts & _targets)
    return docking, n > 0, n


def _score_hl(orig):
    """Mirror code_hl/analyze_replicates.py's scoring block for one compound."""
    from hl_gap_module import scoring_function
    with _quiet_stdout():
        score, _aux = scoring_function(orig)
    return None if (score is None or score >= 100.0) else float(score)


def _metrics(orig):
    """QED/aLogP/SAS/NP via verify_results' parsers (RDKit-only, cheap)."""
    from verify_results import _parse_lipinski, _parse_sas_np
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
    return qed, alogp, sas, np_score


def first_turn_text(sidecar_path, md_path):
    """All proposer text from turn 1: every assistant text segment before the
    first real user prompt after the task prompt (the adversary feedback).

    Source is the JSON sidecar's full message history — resumed sessions only
    log post-resume sections into the .md, so the sidecar is authoritative.
    Handles all three transcript formats:

      molopt_oa.py unified  {'role', 'content': str | [{'type': 'text'|'tool_result', ...}]}
      molopt_oa.py Gemini   {'role': 'user'|'model', 'parts': [{'text': ...} |
                            {'function_call': ...} | {'function_response': ...}]}
      molopt.py (Ollama)    {'role': 'assistant', 'content': str} and
                            {'role': 'tool', 'tool_name': ..., 'content': str}

    Intra-turn machinery that does NOT end the turn: tool results (any format),
    and the two known retry/nudge prompts injected inside chat_turn ("The
    previous call failed with a connection error:" and "You have reached the
    tool-call limit"). The first real user prompt is the task itself (skipped);
    the second ends turn 1. Falls back to cutting the .md at the first
    "# Adversary feedback" section if the sidecar can't be read.
    """
    _TURN_CONTINUES = ("The previous call failed with a connection error:",
                       "You have reached the tool-call limit",
                       # Non-critique adversary-stage messages: the protocol's
                       # empty-response nudge and its adversary-unreachable
                       # fallback. Neither is real adversary feedback, and the
                       # .md logs both under "# Adversary feedback: [...]".
                       "Your last response had no text.",
                       "The adversary model could not be reached")
    try:
        with open(sidecar_path) as f:
            msgs = json.load(f).get('messages') or []
    except Exception:
        try:
            md = open(md_path).read()
        except Exception:
            return ''
        cut = re.search(r'^# Adversary feedback', md, re.M)
        blocks = all_model_response_blocks(md[:cut.start()] if cut else md)
        return '\n'.join(blocks)

    def message_parts(m):
        """(assistant_text, is_tool_result, prompt_text) for one message."""
        role, c, parts = m.get('role'), m.get('content'), m.get('parts') or []
        if role == 'tool':
            return '', True, ''
        texts, tool_result = [], False
        if isinstance(c, str):
            texts.append(c)
        elif isinstance(c, list):
            if any(isinstance(b, dict) and b.get('type') == 'tool_result' for b in c):
                tool_result = True
            texts = [' '.join(b.get('text', '') for b in c
                              if isinstance(b, dict) and b.get('type') == 'text')]
        texts += [p.get('text', '') for p in parts
                  if isinstance(p, dict) and p.get('text')]
        if any(isinstance(p, dict) and 'function_response' in p for p in parts):
            tool_result = True
        return '\n'.join(t for t in texts if t), tool_result, ''

    parts, seen_task_prompt = [], False
    for m in msgs:
        role = m.get('role')
        text, tool_result, _ = message_parts(m)
        if role == 'user':
            if tool_result or (text and text.startswith(_TURN_CONTINUES)):
                continue
            if not seen_task_prompt:
                seen_task_prompt = True  # the task prompt itself
                continue
            break  # second real prompt: adversary feedback
        if role in ('assistant', 'model') and text:
            parts.append(text)
    return '\n'.join(parts)


def load_done(path):
    """(set_label, replicate, canonical_smiles) rows already scored."""
    done = set()
    if os.path.isfile(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                done.add((r['set_label'], r['replicate'], r['canonical_smiles']))
    return done


def append_rows(path, fieldnames, rows):
    exists = os.path.isfile(path)
    with open(path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if not exists:
            w.writeheader()
        w.writerows(rows)


def rewrite_summary(path, stem, study):
    """summary_<stem>.csv and best_per_replicate_<stem>.csv from the compounds CSV."""
    import pandas as pd
    df = pd.read_csv(path)
    metric = 'gap' if study == 'hl' else 'docking'
    cols = [metric, 'qed', 'alogp', 'sas', 'np']
    rows = []
    for label, g in df.groupby('set_label', sort=True):
        rec = {'set_label': label, 'n_replicates': g['replicate'].nunique(),
               'n_compounds': len(g), 'n_unique': g['canonical_smiles'].nunique()}
        if study == 'hl':
            rec['n_gap_failed'] = int(g['gap'].isna().sum())
        else:
            rec['n_dock_failed'] = int(g['docking'].isna().sum())
            rec['n_in_pocket'] = int((g['docked_in_pocket'] == True).sum())  # noqa: E712
        for col in cols:
            s = g[col].dropna()
            rec[f'{col}_mean'] = round(s.mean(), 3) if len(s) else None
            if col == metric:
                rec[f'{col}_median'] = round(s.median(), 3) if len(s) else None
                rec[f'{col}_best'] = round(s.min(), 3) if len(s) else None
        rows.append(rec)
    pd.DataFrame(rows).set_index('set_label').to_csv(path.replace(f'compounds_{stem}', f'summary_{stem}'))

    best = []
    for (label, rep), g in df.groupby(['set_label', 'replicate'], sort=True):
        gd = g.dropna(subset=[metric])
        row = g.iloc[0] if not len(gd) else gd.loc[gd[metric].idxmin()]
        rec = {'set_label': label, 'replicate': rep,
               'canonical_smiles': row['canonical_smiles'], metric: row[metric]}
        if study == 'dock':
            rec['docked_in_pocket'] = row['docked_in_pocket']
            rec['n_target_contacts'] = row['n_target_contacts']
        rec.update({c: row[c] for c in ['qed', 'alogp', 'sas', 'np']})
        best.append(rec)
    pd.DataFrame(best).to_csv(path.replace(f'compounds_{stem}', f'best_per_replicate_{stem}'), index=False)


def main():
    p = argparse.ArgumentParser(description='First-response baseline (one proposer reply, tools allowed).')
    p.add_argument('--study', choices=['dock', 'hl'], required=True)
    args = p.parse_args()

    os.makedirs(OUT_DIR[args.study], exist_ok=True)
    out_path = os.path.join(OUT_DIR[args.study], f'compounds_{STEM}.csv')

    if args.study == 'dock':
        fieldnames = ['set_label', 'replicate', 'proposer_provider', 'proposer_model',
                      'adversary_provider', 'adversary_model', 'protein',
                      'original_smiles', 'canonical_smiles', 'inchikey',
                      'docking', 'docked_in_pocket', 'n_target_contacts',
                      'qed', 'alogp', 'sas', 'np', 'source_turn']
    else:
        fieldnames = ['set_label', 'replicate', 'proposer_provider', 'proposer_model',
                      'adversary_provider', 'adversary_model', 'method',
                      'original_smiles', 'canonical_smiles', 'inchikey',
                      'gap', 'qed', 'alogp', 'sas', 'np', 'source_turn']

    done = load_done(out_path)
    for batch_id, (key, model) in BATCHES[args.study].items():
        manifest = json.load(open(os.path.join(_BATCHES, batch_id, 'manifest.json')))
        provider = 'ollama' if key in ('kimi', 'deepseek') else key
        for e in manifest['entries']:
            rep = str(e['replicate'])
            n_new = 0
            # First proposer response = everything the proposer says in turn 1,
            # which may span several text segments separated by intra-turn tool
            # calls. The sidecar holds the full history even for resumed sessions
            # whose .md lacks the original turn 1; fall back to the .md cut.
            turn1 = first_turn_text(e.get('sidecar_path'), e['md_path'])
            comps = list(extract_smiles(turn1, MIN_HEAVY_ATOMS)) if turn1 else []
            rows = []
            for orig, canon, mol in comps:
                if (key, rep, canon) in done:
                    continue
                from rdkit import Chem
                row = {'set_label': key, 'replicate': rep,
                       'proposer_provider': provider, 'proposer_model': model,
                       'adversary_provider': '', 'adversary_model': '',
                       'original_smiles': orig, 'canonical_smiles': canon,
                       'inchikey': Chem.MolToInchiKey(mol),
                       'qed': None, 'alogp': None, 'sas': None, 'np': None,
                       'source_turn': 0}
                if args.study == 'dock':
                    row['protein'] = e.get('protein', 'HMGCR')
                    try:
                        row['docking'], row['docked_in_pocket'], row['n_target_contacts'] = _score_dock(orig)
                    except Exception as err:
                        print(f'  dock failed {canon}: {err}')
                        row['docking'], row['docked_in_pocket'], row['n_target_contacts'] = None, False, None
                else:
                    row['method'] = 'GFN2-xTB'
                    try:
                        row['gap'] = _score_hl(orig)
                    except Exception as err:
                        print(f'  gap failed {canon}: {err}')
                        row['gap'] = None
                qed, alogp, sas, np_score = _metrics(orig)
                row.update({'qed': qed, 'alogp': alogp, 'sas': sas, 'np': np_score})
                rows.append(row)
            if rows:
                append_rows(out_path, fieldnames, rows)
                n_new = len(rows)
            status = 'no SMILES in first reply' if not comps else f'{n_new} new / {len(comps)} extracted'
            print(f'{batch_id} rep{rep} [{key}]: {status}')
    rewrite_summary(out_path, STEM, args.study)
    print(f'done -> {out_path}')


if __name__ == '__main__':
    main()