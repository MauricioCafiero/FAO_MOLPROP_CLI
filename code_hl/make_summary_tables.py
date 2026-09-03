#!/usr/bin/env python3
"""make_summary_tables.py - Build SUMMARY_TABLES.md for the HL-gap study.

The HL counterpart of results/batches/analysis/SUMMARY_TABLES.md. Same convention:
every metric is computed PER REPLICATE FIRST, THEN AVERAGED across replicates (not
pooled over all compounds), and "Spread" is the SD across replicates.

Metric set is smaller than docking's by design:
  - no in-pocket rows -- there is no pocket for a HOMO-LUMO objective.
  - no QED rows -- QED is a drug-likeness score; it has no meaning for a materials
    target, where the same molecule can be an excellent result and a terrible drug.
  - SAS is kept: synthetic accessibility still matters for something you intend to make.

Both GA pools (frag10 and full) share one table, as they are two configurations of
one non-LLM system rather than two proposers.

Usage: python code_hl/make_summary_tables.py > results/batches/hl_batches/analysis/SUMMARY_TABLES.md
"""
import os
import sys
import glob
import math
import csv
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
HL = os.path.join(_ROOT, 'results', 'batches', 'hl_batches')

# display label -> baseline set_label in the shot CSVs
MAIN5 = [('OpenAI gpt-5.2', 'openai'), ('Anthropic haiku-4.5', 'anthropic'),
         ('Gemini 3-flash', 'gemini'), ('kimi k2.6', 'kimi-k2.6'),
         ('deepseek v4-pro', 'deepseek-v4-pro')]
EXTRA7 = [('gemma4', 'gemma4'), ('glm-5.2', 'glm-5.2'),
          ('nemotron-3-ultra', 'nemotron-3-ultra'), ('nemotron-3-super', 'nemotron-3-super'),
          ('nemotron-3-nano', 'nemotron-3-nano'), ('gpt-oss-20b', 'gpt-oss-20b'),
          ('gpt-oss-120b', 'gpt-oss-120b')]

AGENTIC = [
    ('OpenAI gpt-5.2', 'hl_gpt-5.2_vs_gpt-5.2_5x4'),
    ('Anthropic haiku-4.5', 'hl_claude-haiku-4-5_vs_claude-haiku-4-5_5x4'),
    ('Gemini 3-flash', 'hl_gemini-3-flash-preview_vs_gemini_5x4'),
    ('kimi k2.6', 'hl_kimi-k2.6_vs_kimi-k2.6_5x4'),
    ('deepseek v4-pro', 'hl_deepseek-v4-pro_vs_deepseek-v4-pro_5x4'),
]
GA = [('GA restricted (frag10)', '5x4'), ('GA unrestricted (full pool)', '5x4_full')]


def _rows(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def _f(row, key):
    try:
        return float(row[key])
    except (TypeError, ValueError, KeyError):
        return None


def per_rep(rows, label_filter=None):
    """{replicate: [rows]} for one model, dropping failed-gap rows."""
    out = defaultdict(list)
    for r in rows:
        if label_filter is not None and r.get('set_label') != label_filter:
            continue
        if _f(r, 'gap') is None:
            continue
        out[r['replicate']].append(r)
    return out


def mean_sd(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None, None
    m = sum(vals) / len(vals)
    if len(vals) < 2:
        return m, 0.0
    var = sum((v - m) ** 2 for v in vals) / (len(vals) - 1)
    return m, math.sqrt(var)


def fmt(v, nd=3):
    return '—' if v is None else f'{v:.{nd}f}'


def stats_for(reps, tool_rows=None):
    tops, means, sas = [], [], []
    for rep, rs in reps.items():
        gaps = [_f(r, 'gap') for r in rs]
        gaps = [g for g in gaps if g is not None]
        if not gaps:
            continue
        tops.append(min(gaps))
        means.append(sum(gaps) / len(gaps))
        s = [_f(r, 'sas') for r in rs]
        s = [x for x in s if x is not None]
        if s:
            sas.append(sum(s) / len(s))
    t_m, t_sd = mean_sd(tops)
    g_m, g_sd = mean_sd(means)
    s_m, s_sd = mean_sd(sas)
    turns = calls = None
    if tool_rows is not None:
        tv = [_f(r, 'n_turns') for r in tool_rows]
        cv = [_f(r, 'total_tool_calls') for r in tool_rows]
        turns, _ = mean_sd(tv)
        calls, _ = mean_sd(cv)
    return {'n': len(tops), 'top': t_m, 'top_sd': t_sd, 'gap': g_m, 'gap_sd': g_sd,
            'sas': s_m, 'sas_sd': s_sd, 'turns': turns, 'calls': calls}


def table(headers, cols, static_turns=None, static_calls=None):
    """cols: list of stat dicts aligned with headers."""
    lines = ['| Metric | ' + ' | '.join(headers) + ' |',
             '|---' + '|---:' * len(headers) + '|']
    def row(name, key, nd=3):
        return f'| {name} | ' + ' | '.join(fmt(c[key], nd) for c in cols) + ' |'
    lines.append(row('Average top gap (eV)', 'top'))
    lines.append(row('Spread of top gap (SD)', 'top_sd'))
    lines.append(row('Average of average gap (eV)', 'gap'))
    lines.append(row('Spread of average gap (SD)', 'gap_sd'))
    lines.append(row('Average of average SAS', 'sas', 2))
    lines.append(row('Spread of average SAS (SD)', 'sas_sd', 2))
    if static_turns is not None:
        lines.append('| Average number of turns | ' + ' | '.join([static_turns] * len(cols)) + ' |')
        lines.append('| Average number of tool calls | ' + ' | '.join([static_calls] * len(cols)) + ' |')
    else:
        lines.append(row('Average number of turns', 'turns', 1))
        lines.append(row('Average number of tool calls', 'calls', 1))
    lines.append('')
    lines.append('n (replicates): ' + ', '.join(f'{h} {c["n"]}' for h, c in zip(headers, cols)))
    return '\n'.join(lines)


def main():
    print('# Summary Tables — HOMO-LUMO Gap Study (GFN2-xTB)\n')
    print('Per-model summary tables, one section per protocol. Each metric is computed '
          '**per replicate first, then averaged across replicates** (not pooled across all '
          'compounds) — e.g. "average of average gap" is the mean of each replicate\'s own '
          'mean gap, not the mean over every compound regardless of replicate. "Spread" = '
          'standard deviation across replicates. Lower gap is better throughout.\n')
    print('This table set is deliberately smaller than the HMGCR study\'s. There is no '
          'in-pocket row (a HOMO-LUMO objective has no binding pocket), and no QED row '
          '(QED scores drug-likeness, which carries no meaning for a materials target — a '
          'molecule can be an excellent small-gap result and a hopeless drug). SAS is kept: '
          'synthetic accessibility still matters for anything you intend to make. Gap and '
          'SAS are therefore the two chemistry metrics reported.\n')
    print('Compounds whose GFN2-xTB gap calculation failed are excluded before averaging, '
          'so a failed calculation is never counted as a real (large) gap.\n')

    for shot, title, desc in (
        ('zero', 'Zero-Shot',
         'single API call per replicate, no tools, no adversary, no multi-turn refinement. '
         'The user message names only the objective (`GFN2-xTB HOMO-LUMO gap`).'),
        ('frag', 'Fragment-Suggested Zero-Shot',
         'as zero-shot, but the system prompt also lists the 9 base rings and the exact 10 '
         'fragments `adversarial_set.md` was enumerated from.'),
        ('few', 'Few-Shot',
         'single API call per replicate; the user message is the full 260-molecule '
         '`adversarial_set.md` gap list, with "learn the trends, then propose" instructions.'),
    ):
        rows = _rows(os.path.join(HL, f'{shot}_shot', 'analysis', f'compounds_{shot}_shot.csv'))
        print(f'## {title}\n')
        print(f'**Protocol:** {desc} 5 replicates/model.\n')
        cols = [stats_for(per_rep(rows, lbl)) for _, lbl in MAIN5]
        print(table([h for h, _ in MAIN5], cols, static_turns='1.0', static_calls='0.0'))
        print()
        print(f'### {title} — additional Ollama models\n')
        print('Same protocol, for the 7 Ollama-cloud-only models that have no agentic-loop '
              'role in this study (baselines only).\n')
        cols7 = [stats_for(per_rep(rows, lbl)) for _, lbl in EXTRA7]
        print(table([h for h, _ in EXTRA7], cols7, static_turns='1.0', static_calls='0.0'))
        print()

    print('## Agentic 5×4 (self-critique)\n')
    print('**Protocol:** proposer designs, calls chemistry/gap tools, and is critiqued by '
          'itself (same model as adversary) for up to 5 refinement turns with a 4 tool-round '
          'cap per turn. 5 replicates/model.\n')
    cols, heads = [], []
    for label, batch in AGENTIC:
        comp = _rows(os.path.join(HL, batch, 'analysis', f'compounds_{batch}.csv'))
        tool = _rows(os.path.join(HL, batch, 'analysis', f'tool_usage_{batch}.csv'))
        cols.append(stats_for(per_rep(comp), tool_rows=tool))
        heads.append(label)
    print(table(heads, cols))
    print()
    print('A gemma4 self-critique set also exists under `results/batches/hl_batches/` but was '
          'a test run, not a study condition, and is excluded here and from the statistics.\n')

    print('## Genetic Algorithm baselines (non-LLM)\n')
    print('**Protocol:** population-based search over the same 9 base rings, fitness = the '
          'real GFN2-xTB gap, no LLM and no chemical reasoning. pop=5 × 4 generations, 5 '
          'replicates each. Both pools share one table: they are two configurations of one '
          'system, not two proposers. `frag10` is restricted to the exact 10 fragments '
          'frag-shot shows the model (the pool-matched comparison for frag-shot); `full` '
          'searches the entire ~390-item substituent pool (the comparison for zero-shot, '
          'which is shown no fragment menu).\n')
    gcols, gheads = [], []
    for label, sub in GA:
        p = os.path.join(HL, 'ga_baseline', sub, 'analysis', f'compounds_ga_{sub}.csv')
        gcols.append(stats_for(per_rep(_rows(p))))
        gheads.append(label)
    print(table(gheads, gcols, static_turns='—', static_calls='—'))
    print()
    print('Turns and tool calls are not applicable: the GA has no conversation and calls no '
          'tools — it evaluates genomes directly.\n')


if __name__ == '__main__':
    sys.exit(main())
