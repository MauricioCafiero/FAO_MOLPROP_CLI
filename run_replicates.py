#!/usr/bin/env python3
"""
run_replicates.py - Run N replicates of each adversary set, sequentially.

Adversary sets (each is one (proposer, adversary) pairing):
  openai_vs_anthropic   molopt_oa.py --start openai     (OpenAI proposes+tools, Anthropic critiques)
  openai_vs_openai      molopt_oa.py --start openai --adversary openai (OpenAI proposes+tools,
                          OpenAI self-critiques -- same model both sides; matches the gpt-5.2
                          adversary used by the ollama_*/gemini_* sets for cross-set comparability)
  anthropic_vs_openai   molopt_oa.py --start anthropic  (Anthropic proposes+tools, OpenAI critiques)
  ollama_vs_openai      molopt.py   --adversary openai   (Ollama proposes+tools, OpenAI critiques)
  ollama_vs_anthropic   molopt.py   --adversary anthropic (Ollama proposes+tools, Anthropic critiques)
  gemini_vs_openai      molopt_oa.py --start gemini --adversary openai     (Gemini proposes+tools, OpenAI critiques)
  gemini_vs_anthropic   molopt_oa.py --start gemini --adversary anthropic  (Gemini proposes+tools, Anthropic critiques)

For each set, runs --replicates sessions back-to-back. Docking is CPU-bound so runs are
strictly sequential (no parallelism). Each replicate gets its own --results-dir under
results/batches/<batch_id>/<set>/rep<N>/, so one run -> one .md + one .json sidecar.

A manifest (manifest.json at the batch root) records every job; analyze_replicates.py reads it.

Resumable: a replicate with a terminal sidecar (status Done / max_turns_reached) is skipped on
re-launch, so a killed batch loses no completed work. --force re-runs everything. We never
delete run artifacts -- a re-run just writes a new timestamped .md; the manifest points at the
latest terminal one.

Env: the venv MUST be activated (source fao-env/bin/activate) so fao-env/bin is on PATH -- the
openbabel Python bindings and the obabel CLI (used by dockstring) both need it. This script
refuses to run otherwise. API keys are inherited from the environment (source ~/.zshrc first).

Launch detached so a long batch survives the shell (per the long-job pattern in CLAUDE.md):
  source fao-env/bin/activate
  fao-env/bin/python -c "import subprocess; subprocess.Popen(['fao-env/bin/python','run_replicates.py','--replicates','3'], stdout=open('run_replicates.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)"
Monitor: tail -f run_replicates.log ; or read manifest.json.
"""

import os
import sys
import json
import glob
import shutil
import time
import argparse
import subprocess

_HERE = os.path.dirname(os.path.abspath(__file__))
_PY = sys.executable  # same interpreter (the activated venv's python)
_DEFAULT_RESULTS_ROOT = os.path.join(_HERE, 'results', 'batches')

# Adversary set definitions: label -> (script, argv template). Model names are substituted
# from the CLI flags at build time so all sets share one --openai-model / --anthropic-model /
# --ollama-model.
def build_configs(models):
    """Return an ordered dict label -> {script, args} using the resolved model names."""
    om, oai, ant, gem = models['ollama'], models['openai'], models['anthropic'], models['gemini']
    return {
        'openai_vs_anthropic': {
            'script': 'molopt_oa.py',
            'args': ['--start', 'openai', '--openai-model', oai, '--anthropic-model', ant],
        },
        'openai_vs_openai': {
            'script': 'molopt_oa.py',
            'args': ['--start', 'openai', '--adversary', 'openai', '--openai-model', oai],
        },
        'anthropic_vs_openai': {
            'script': 'molopt_oa.py',
            'args': ['--start', 'anthropic', '--openai-model', oai, '--anthropic-model', ant],
        },
        'ollama_vs_openai': {
            'script': 'molopt.py',
            'args': ['--model', om, '--adversary', 'openai', '--adversary-model', oai],
        },
        'ollama_vs_anthropic': {
            'script': 'molopt.py',
            'args': ['--model', om, '--adversary', 'anthropic', '--adversary-model', ant],
        },
        'gemini_vs_openai': {
            'script': 'molopt_oa.py',
            'args': ['--start', 'gemini', '--adversary', 'openai',
                     '--gemini-model', gem, '--openai-model', oai],
        },
        'gemini_vs_anthropic': {
            'script': 'molopt_oa.py',
            'args': ['--start', 'gemini', '--adversary', 'anthropic',
                     '--gemini-model', gem, '--anthropic-model', ant],
        },
    }


# Sidecar statuses that mean a run finished (proposer said Done or hit the turn cap).
_TERMINAL = {'Done', 'max_turns_reached'}


def sidecar_status(repdir):
    """Return (status, sidecar_path, md_path) for the terminal sidecar in repdir, else (None, None, None).

    A rep is 'complete' if ANY sidecar in its dir has a terminal status. The matching .md is the
    sidecar's sibling (same stem). Multiple .md can accumulate from re-runs; we want the terminal
    one, which is the one analyze_replicates.py should read.
    """
    jsons = sorted(glob.glob(os.path.join(repdir, '*.json')),
                   key=os.path.getmtime, reverse=True)
    for jp in jsons:
        try:
            with open(jp) as f:
                payload = json.load(f)
        except Exception:
            continue
        if payload.get('status') in _TERMINAL:
            md = os.path.splitext(jp)[0] + '.md'
            return payload['status'], jp, md
    return None, None, None


def newest_run_files(repdir):
    """After a fresh run, return (sidecar_path, md_path) for the newest sidecar in repdir."""
    jsons = sorted(glob.glob(os.path.join(repdir, '*.json')),
                   key=os.path.getmtime, reverse=True)
    if not jsons:
        return None, None
    jp = jsons[0]
    return jp, os.path.splitext(jp)[0] + '.md'


def load_manifest(path):
    if os.path.isfile(path):
        try:
            with open(path) as f:
                return json.load(f)
        except Exception:
            pass
    return {'entries': []}


def write_manifest(path, manifest):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp, path)


def build_argv(cfg, *, repdir, protein, max_turns, max_tool_calls, quiet, mock_tools):
    """Full argv for one replicate run."""
    argv = [_PY, os.path.join(_HERE, cfg['script'])]
    argv += list(cfg['args'])
    argv += ['--protein', protein,
             '--results-dir', repdir,
             '--max-turns', str(max_turns),
             '--max-tool-calls', str(max_tool_calls)]
    if quiet:
        argv.append('--quiet')
    if mock_tools:
        argv.append('--mock-tools')
    return argv


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog='run_replicates.py',
        description='Run N replicates of each adversary set sequentially; writes a manifest '
                    'for analyze_replicates.py.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--replicates', type=int, default=3,
                   help='Replicates per adversary set (default: 3).')
    p.add_argument('--sets',
                   default='openai_vs_anthropic,anthropic_vs_openai,ollama_vs_openai,ollama_vs_anthropic',
                   help='Comma-separated subset of adversary sets to run (default: the 4 '
                        'non-Gemini sets). Choices: openai_vs_anthropic, openai_vs_openai, '
                        'anthropic_vs_openai, ollama_vs_openai, ollama_vs_anthropic, '
                        'gemini_vs_openai, gemini_vs_anthropic.')
    p.add_argument('--protein', default='HMGCR', help='Docking target (default: HMGCR).')
    p.add_argument('--max-turns', type=int, default=20,
                   help='--max-turns passed to each run (default: 20).')
    p.add_argument('--max-tool-calls', type=int, default=12,
                   help='--max-tool-calls passed to each run (default: 12).')
    p.add_argument('--ollama-model', default='deepseek-v4-pro',
                   help='Ollama main model for the ollama_* sets (default: deepseek-v4-pro).')
    p.add_argument('--openai-model', default='gpt-5.2',
                   help='OpenAI model (default: gpt-5.2).')
    p.add_argument('--anthropic-model', default='claude-haiku-4-5-20251001',
                   help='Anthropic model (default: claude-haiku-4-5-20251001).')
    p.add_argument('--gemini-model', default='gemini-3-flash-preview',
                   help='Gemini model for the gemini_* sets (default: gemini-3-flash-preview).')
    p.add_argument('--batch-id', default=None,
                   help='Batch directory name (default: <protein>_<timestamp>). '
                        'All output goes under results/batches/<batch-id>/.')
    p.add_argument('--results-root', default=_DEFAULT_RESULTS_ROOT,
                   help='Where batch dirs are written (default: results/batches).')
    p.add_argument('--force', action='store_true',
                   help='Re-run every replicate, ignoring completed ones (does not delete old '
                        'artifacts; the manifest is updated to the newest run).')
    p.add_argument('--verbose-runs', action='store_true',
                   help='Do NOT pass --quiet to the runs (default: runs are quiet; the .md is '
                        'the record).')
    p.add_argument('--dry-run', action='store_true',
                   help='Print the command for each replicate and exit without running.')
    p.add_argument('--mock-tools', action='store_true',
                   help='Pass --mock-tools through to every run: replaces the docking-dependent '
                        'tools with an instant synthetic score, skipping real Vina docking. For '
                        'smoke-testing a set/wiring change in seconds-per-rep instead of '
                        'minutes/hours -- not for real runs.')
    args = p.parse_args(argv)

    # --- Env guard: the whole stack only works with fao-env/bin on PATH. ---
    if shutil.which('obabel') is None:
        print("ERROR: 'obabel' not found on PATH. Activate the venv first:\n"
              "  source fao-env/bin/activate\n"
              "(openbabel bindings + the obabel CLI dockstring shells out to both need\n"
              " fao-env/bin on PATH; fao-env/bin/python alone is not enough.)",
              file=sys.stderr)
        return 2

    configs = build_configs({'ollama': args.ollama_model, 'openai': args.openai_model,
                             'anthropic': args.anthropic_model, 'gemini': args.gemini_model})
    selected = [s.strip() for s in args.sets.split(',') if s.strip()]
    bad = [s for s in selected if s not in configs]
    if bad:
        print(f"ERROR: unknown set(s): {bad}. Choices: {list(configs)}", file=sys.stderr)
        return 2

    batch_id = args.batch_id or f"{args.protein}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
    batch_dir = os.path.join(args.results_root, batch_id)
    manifest_path = os.path.join(batch_dir, 'manifest.json')
    os.makedirs(batch_dir, exist_ok=True)

    manifest = load_manifest(manifest_path)
    manifest.setdefault('entries', [])
    manifest['batch_id'] = batch_id
    manifest['batch_dir'] = batch_dir
    manifest['protein'] = args.protein
    manifest['replicates'] = args.replicates
    manifest['sets'] = selected
    manifest['mock_tools'] = args.mock_tools
    manifest['models'] = {'ollama': args.ollama_model, 'openai': args.openai_model,
                          'anthropic': args.anthropic_model, 'gemini': args.gemini_model}
    write_manifest(manifest_path, manifest)

    total = len(selected) * args.replicates
    done = 0
    print(f"Batch {batch_id}: {len(selected)} set(s) x {args.replicates} reps = {total} runs "
          f"(sequential). Manifest -> {manifest_path}")
    if args.dry_run:
        print("(dry-run: commands only)\n")

    for label in selected:
        cfg = configs[label]
        set_dir = os.path.join(batch_dir, label)
        os.makedirs(set_dir, exist_ok=True)
        for rep in range(1, args.replicates + 1):
            done += 1
            repdir = os.path.join(set_dir, f'rep{rep}')
            os.makedirs(repdir, exist_ok=True)
            argv = build_argv(cfg, repdir=repdir, protein=args.protein,
                              max_turns=args.max_turns, max_tool_calls=args.max_tool_calls,
                              quiet=not args.verbose_runs, mock_tools=args.mock_tools)

            # Skip if already complete (unless --force).
            if not args.force:
                st, jp, md = sidecar_status(repdir)
                if st is not None:
                    print(f"[{done}/{total}] {label} rep{rep}: SKIP (already {st})")
                    _upsert(manifest, batch_id, label, rep, cfg['script'], args.protein,
                            'skipped', md_path=md, sidecar_path=jp, exit_code=None)
                    write_manifest(manifest_path, manifest)
                    continue

            if args.dry_run:
                print(f"[{done}/{total}] {label} rep{rep} -> {' '.join(argv)}")
                continue

            print(f"[{done}/{total}] {label} rep{rep}: RUN ...")
            _upsert(manifest, batch_id, label, rep, cfg['script'], args.protein,
                    'running', md_path=None, sidecar_path=None, exit_code=None)
            write_manifest(manifest_path, manifest)
            t0 = time.time()
            try:
                proc = subprocess.run(argv, cwd=_HERE)
                exit_code = proc.returncode
            except Exception as err:
                print(f"  [launch failed: {err}]", file=sys.stderr)
                exit_code = None
            elapsed = time.time() - t0

            jp, md = newest_run_files(repdir)
            st = None
            if jp:
                try:
                    with open(jp) as f:
                        st = json.load(f).get('status')
                except Exception:
                    st = None
            if st in _TERMINAL and exit_code == 0:
                status = 'complete'
            elif exit_code == 0 and st is None:
                status = 'incomplete'  # ran but no sidecar (e.g. crashed before first write)
            else:
                status = 'failed'
            print(f"  -> {status} (exit={exit_code}, sidecar={st}, {elapsed:.0f}s)")
            _upsert(manifest, batch_id, label, rep, cfg['script'], args.protein,
                    status, md_path=md, sidecar_path=jp, exit_code=exit_code)
            write_manifest(manifest_path, manifest)

    write_manifest(manifest_path, manifest)
    n_complete = sum(1 for e in manifest['entries'] if e['status'] == 'complete')
    n_skipped = sum(1 for e in manifest['entries'] if e['status'] == 'skipped')
    n_failed = sum(1 for e in manifest['entries'] if e['status'] == 'failed')
    print(f"\nBatch {batch_id} done. complete={n_complete} skipped={n_skipped} "
          f"failed={n_failed}. Manifest: {manifest_path}")
    print(f"Analyze with: fao-env/bin/python analyze_replicates.py --batch-dir {batch_dir}")
    return 0


def _upsert(manifest, batch_id, label, rep, script, protein, status,
            md_path, sidecar_path, exit_code):
    """Insert or update the manifest entry for (label, rep)."""
    entries = manifest['entries']
    for e in entries:
        if e['set_label'] == label and e['replicate'] == rep:
            entry = e
            break
    else:
        entry = {'set_label': label, 'replicate': rep}
        entries.append(entry)
    entry.update({'batch_id': batch_id, 'script': script, 'protein': protein,
                  'status': status, 'md_path': md_path, 'sidecar_path': sidecar_path,
                  'exit_code': exit_code})


if __name__ == '__main__':
    sys.exit(main())