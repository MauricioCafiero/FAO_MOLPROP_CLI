#!/usr/bin/env python3
"""
run_zero_few_shot.py - Non-agentic baseline for the HOMO-LUMO gap objective:
single-shot molecule generation (no tools, no adversary, no multi-turn
refinement) across the same proposers used elsewhere in this study, for
zero-shot, few-shot, and fragment-suggested zero-shot prompting.

The HL counterpart of code/run_zero_few_shot.py (the HMGCR docking version).
Structure, CLI, output layout and provider callers are ported unchanged; only
the task changes, so the two objectives stay directly comparable:

  zero-shot: system prompt only states the task; the user message names just
             the objective ('GFN2-xTB HOMO-LUMO gap'), the direct analogue of
             the docking version's bare protein name.
  few-shot:  system prompt is hl_gap_module.task_specific_prompt (the same text
             the agentic HL runs use) plus the "learn trends, then propose"
             scaffold; the user message is the full code_hl/adversarial_set.md
             molecule/gap list, worded exactly as molopt.py's first_prompt so
             the few-shot condition sees what turn 1 of an agentic run sees.
  frag-shot: same task framing and same bare-objective user message as
             zero-shot, but the system prompt also suggests base rings and
             functional groups to build from.

The frag-shot functional groups are the exact 10 substituents that
code_hl/adversarial_set.md was enumerated from (sub_cycle over base_rings x
clean_ring_locations x these 10 = 26 x 10 = 260 entries), so frag-shot, the
seed set, and a pool-matched GA baseline all share one fragment vocabulary.
Note these are all linker+group combinations, e.g. 'N(I)' is -NH- + iodo, not
a bare iodo -- the docking version's menu differs here (bare I / C#N).

Each (shot-type, model, replicate) is a single API call. Output matches the
agentic HL convention (method GFN2-xTB, <model>_GFN2-xTB_<ts>.md + .json
sidecar, manifest.json at the <shot>_shot level) so one analyzer can read
baselines and agentic batches alike.

Usage:
  python3 code_hl/run_zero_few_shot.py --shot zero --replicates 5
  python3 code_hl/run_zero_few_shot.py --shot few  --replicates 5
  python3 code_hl/run_zero_few_shot.py --shot frag --replicates 5
  python3 code_hl/run_zero_few_shot.py --shot zero --replicates 1 --models openai,gemini
"""
import os
import sys
import json
import time
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))  # code_hl/
_ROOT = os.path.dirname(_HERE)  # repo root
sys.path.insert(0, _HERE)

from hl_gap_module import task_specific_prompt  # noqa: E402  (the HL few-shot prompt)

_DEFAULT_API_TIMEOUT = 120.0
METHOD = 'GFN2-xTB'

OPENAI_DEFAULT_MODEL = 'gpt-5.2'
ANTHROPIC_DEFAULT_MODEL = 'claude-haiku-4-5-20251001'
GEMINI_DEFAULT_MODEL = 'gemini-3-flash-preview'
OLLAMA_HOST = 'https://ollama.com'

# label -> (provider, default_model). Identical roster to the docking baseline.
MODELS = {
    'openai': ('openai', OPENAI_DEFAULT_MODEL),
    'anthropic': ('anthropic', ANTHROPIC_DEFAULT_MODEL),
    'gemini': ('gemini', GEMINI_DEFAULT_MODEL),
    'kimi-k2.6': ('ollama', 'kimi-k2.6'),
    'deepseek-v4-pro': ('ollama', 'deepseek-v4-pro'),
    'gemma4': ('ollama', 'gemma4'),
    'glm-5.2': ('ollama', 'glm-5.2'),
    'nemotron-3-ultra': ('ollama', 'nemotron-3-ultra'),
    'nemotron-3-super': ('ollama', 'nemotron-3-super'),
    'nemotron-3-nano': ('ollama', 'nemotron-3-nano:30b'),
    'gpt-oss-20b': ('ollama', 'gpt-oss:20b'),
    'gpt-oss-120b': ('ollama', 'gpt-oss:120b'),
}

# --- Prompts -----------------------------------------------------------------
# Structural port of the docking prompts: "best possible docking score (the most
# negative) to a particular protein target" -> "best possible HOMO-LUMO gap (the
# smallest) computed with a particular method". Direction is stated explicitly
# because the agentic HL arm minimises the gap; the objective must match it.
ZERO_SHOT_SYSTEM = '''# You are a materials science assistant. Your task is to design new molecules
with the best possible HOMO-LUMO gap (the smallest), computed with a particular method, given in the first user message.
You will deliver up to five potential molecules in SMILES format, along with reasoning for why you chose those molecules
and an estimate of their HOMO-LUMO gaps.
'''
ZERO_SHOT_USER = 'GFN2-xTB HOMO-LUMO gap'

# zero-shot, but with a fixed set of suggested rings/functional groups to build
# from -- same task framing as ZERO_SHOT_SYSTEM, same bare-objective user
# message, just with fragment suggestions appended to the system prompt.
# Rings are base_rings from MolPropOp.py; groups are the exact 10 substituents
# code_hl/adversarial_set.md was enumerated from.
FRAG_SHOT_SYSTEM = '''# You are a materials science assistant. Your task is to design new molecules
with the best possible HOMO-LUMO gap (the smallest), computed with a particular method, given in the first user message.
You will deliver up to five potential molecules in SMILES format, along with reasoning for why you chose those molecules
and an estimate of their HOMO-LUMO gaps.

## The following are SMILES for rings that you should use as the base of your molecules:
- 'c1ccccc1', #benzene
- 'n1ccccc1', #pyridine
- 'o1cccc1',  #furan
- 's1cccc1',  #thiophene
- '[nH]1cccc1', #pyrrole
- 'n1c[nH]cc1', #imidazole
- 'c1ccc2ccccc2c1', #naphthalene
- 'c1ccc2cc3ccccc3cc2c1', #anthracene
- 'O=c1cc(-c2ccccc2)oc2ccccc12' #flavone

## The following are SMILES for functional groups that you may use to modify the rings; you may also choose to use other functional groups:
- N(I)
- O(C#N)
- C(=O)O(C(C)C)
- C#C(SC)
- C(C(=O)[O-])
- C(C)
- C=C([N+](=O)[O-])
- C(N)
- C([O-])
- CC(N(C)C)
'''
FRAG_SHOT_USER = 'GFN2-xTB HOMO-LUMO gap'

# Same scaffold as the docking few-shot prompt, with hl_gap_module's task text
# (which already states "smallest possible HOMO-LUMO gap") in place of the
# docking one, and "score" -> "gap" throughout.
FEW_SHOT_SYSTEM = f'''
{task_specific_prompt}

## You will first:
- Read the list of molecule SMILES and gaps
- Ascertain any features of the molecules that contribute to a smaller gap. For example, if,
from one molecule to the next, the addition of an O group makes the gap smaller.
- Gather all of these trends across all of the molecules.

## Once you have ascertained the trends:
- Use the trends you learned to suggest 1-5 new molecules that obey the trends you found
and which should have a smaller gap than the molecules in the list.
- Provide reasoning as to why you created those new molecules.
- Estimate the new gaps.
'''


def _few_shot_user() -> str:
    # Wording matches molopt.py's first_prompt exactly, so the few-shot condition
    # receives the identical framing turn 1 of an agentic run receives.
    with open(os.path.join(_HERE, 'adversarial_set.md'), 'r') as f:
        context = f.read()
    return f'\n  Here is a list of molecules and their HOMO-LUMO gaps:\n  {context}\n'


# --- Minimal .env loader (verbatim pattern from molopt.py / molopt_oa.py) ---

_ENV_KEYS = ('OLLAMA_API_KEY', 'OPENAI_API_KEY', 'ANTHROPIC_API_KEY', 'GEMINI_API_KEY')


def load_dotenv(path: str = '.env') -> None:
    if not os.path.isfile(path):
        return
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            key, _, val = line.partition('=')
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            if key in _ENV_KEYS and key not in os.environ:
                os.environ[key] = val


# --- Single-shot callers, one per provider (no tools, no loop) --------------
# Ported unchanged from the docking baseline: lazy import, max_retries=0 +
# timeout so a permanent 401/429 fails fast instead of hanging on SDK backoff.

def call_openai(model, system, user, api_key, timeout=_DEFAULT_API_TIMEOUT):
    from openai import OpenAI
    client = OpenAI(api_key=api_key, timeout=timeout, max_retries=0)
    resp = client.responses.create(model=model, instructions=system, input=user)
    return resp.output_text or ''


def call_anthropic(model, system, user, api_key, timeout=_DEFAULT_API_TIMEOUT,
                    base_url="https://api.anthropic.com"):
    from anthropic import Anthropic
    client = Anthropic(api_key=api_key, base_url=base_url, timeout=timeout, max_retries=0)
    resp = client.messages.create(
        model=model, max_tokens=4096, system=system,
        messages=[{"role": "user", "content": user}])
    return "".join(b.text for b in resp.content if getattr(b, "type", None) == "text").strip()


def call_gemini(model, system, user, api_key, timeout=_DEFAULT_API_TIMEOUT):
    from google import genai
    from google.genai import types
    client = genai.Client(api_key=api_key,
                          http_options=types.HttpOptions(timeout=float(timeout) * 1000))
    cfg = types.GenerateContentConfig(system_instruction=system)
    r = client.models.generate_content(
        model=model,
        contents=[types.Content(role='user', parts=[types.Part(text=user)])],
        config=cfg)
    return ''.join(getattr(p, 'text', '') or '' for p in r.candidates[0].content.parts).strip()


# think=True is impractical for a single-shot call on these labels (kimi-k2.6
# shows 2-23min+ unpredictable latency; the Ollama-only labels are
# uncharacterised). Same set as the docking baseline.
NO_THINK_LABELS = {
    'kimi-k2.6', 'gemma4', 'glm-5.2', 'nemotron-3-ultra', 'nemotron-3-super',
    'nemotron-3-nano', 'gpt-oss-20b', 'gpt-oss-120b',
}


def call_ollama(model, system, user, host, headers, think=True,
                 timeout=_DEFAULT_API_TIMEOUT):
    from ollama import Client as OllamaClient
    client = OllamaClient(host=host, headers=headers, timeout=timeout)
    resp = client.chat(
        model=model,
        messages=[{"role": "system", "content": system}, {"role": "user", "content": user}],
        think=think,
    )
    return resp.message.content or ''


def call_model(label, system, user, keys) -> str:
    provider, model = MODELS[label]
    if provider == 'openai':
        return call_openai(model, system, user, keys['openai'])
    if provider == 'anthropic':
        return call_anthropic(model, system, user, keys['anthropic'])
    if provider == 'gemini':
        return call_gemini(model, system, user, keys['gemini'])
    if provider == 'ollama':
        headers = {'Authorization': f'Bearer {keys["ollama"]}'} if keys['ollama'] else {}
        think = label not in NO_THINK_LABELS
        return call_ollama(model, system, user, OLLAMA_HOST, headers, think=think)
    raise ValueError(f'unknown provider for {label!r}: {provider!r}')


# --- Output: .md + JSON sidecar matching the agentic HL format ---------------
# Filename and 'method' key mirror molopt.py/molopt_oa.py's HL output so a
# single analyzer reads both; the '# Initial model response:' section is what
# the SMILES extractor looks for.

def write_run(results_dir, label, model, system, user, response, shot):
    os.makedirs(results_dir, exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    safe_model = model.replace(':', '').replace('/', '-')
    md_path = os.path.join(results_dir, f"{safe_model}_{METHOD}_{timestamp}.md")
    sidecar_path = os.path.splitext(md_path)[0] + '.json'

    with open(md_path, 'w') as f:
        f.write(f'# {shot.capitalize()}-Shot Design Session - {timestamp}\n')
        f.write(f'# method: {METHOD} | main model: {model} (think=n/a) | adversary: none\n\n')
        f.write('# Initial model response:\n')
        f.write(response + '\n\n')
        f.write('# Session end: Done\n')

    messages = [
        {'role': 'system', 'content': system},
        {'role': 'user', 'content': user},
        {'role': 'assistant', 'content': response},
    ]
    payload = {
        'method': METHOD, 'shot': shot, 'label': label, 'model': model,
        'written_at_turn': 0, 'status': 'Done', 'messages': messages,
    }
    with open(sidecar_path, 'w') as f:
        json.dump(payload, f, indent=2, default=str)

    return md_path, sidecar_path


# --- Batch runner -------------------------------------------------------------

def run_batch(shot, replicates, labels, results_root, quiet=False):
    load_dotenv()
    keys = {
        'openai': os.environ.get('OPENAI_API_KEY') or '',
        'anthropic': os.environ.get('ANTHROPIC_API_KEY') or '',
        'gemini': os.environ.get('GEMINI_API_KEY') or '',
        'ollama': os.environ.get('OLLAMA_API_KEY') or os.environ.get('OLLAMA_KEY') or '',
    }
    if shot == 'zero':
        system, user = ZERO_SHOT_SYSTEM, ZERO_SHOT_USER
    elif shot == 'frag':
        system, user = FRAG_SHOT_SYSTEM, FRAG_SHOT_USER
    else:
        system, user = FEW_SHOT_SYSTEM, _few_shot_user()

    batch_dir = os.path.join(results_root, f'{shot}_shot')
    os.makedirs(batch_dir, exist_ok=True)
    manifest_path = os.path.join(batch_dir, 'manifest.json')
    manifest = {'entries': []}
    if os.path.isfile(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
    done = {(e['set_label'], e['replicate']) for e in manifest['entries']
            if e.get('status') == 'complete'}

    for label in labels:
        provider, model = MODELS[label]
        set_dir = os.path.join(batch_dir, label)
        for rep in range(1, replicates + 1):
            if (label, rep) in done:
                print(f"[{shot}-shot] {label} rep{rep}: SKIP (already complete)")
                continue
            repdir = os.path.join(set_dir, f'rep{rep}')
            print(f"[{shot}-shot] {label} rep{rep}: calling {provider}/{model} ...")
            try:
                response = call_model(label, system, user, keys)
                status, exit_code = 'complete', 0
            except Exception as err:
                print(f"  [error: {err}]")
                response = f"ERROR: {err}"
                status, exit_code = 'failed', 1
            md_path, sidecar_path = write_run(repdir, label, model, system, user, response, shot)
            manifest['entries'].append({
                'set_label': label, 'replicate': rep, 'batch_id': f'{shot}_shot',
                'script': 'run_zero_few_shot.py', 'method': METHOD,
                'status': status, 'md_path': os.path.abspath(md_path),
                'sidecar_path': os.path.abspath(sidecar_path), 'exit_code': exit_code,
            })
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)

    manifest['batch_id'] = f'{shot}_shot'
    manifest['batch_dir'] = os.path.abspath(batch_dir)
    manifest['method'] = METHOD
    manifest['replicates'] = replicates
    manifest['sets'] = labels
    manifest['models'] = {label: MODELS[label][1] for label in labels}
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\n{shot}-shot batch done. Manifest: {manifest_path}")
    print(f"Analyze with: fao-env/bin/python code_hl/analyze_replicates.py --batch-dir {batch_dir}")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog='run_zero_few_shot.py',
        description='Non-agentic zero/few/frag-shot HOMO-LUMO gap baseline.')
    p.add_argument('--shot', choices=['zero', 'few', 'frag'], required=True,
                   help='Which prompting mode to run.')
    p.add_argument('--replicates', type=int, default=5,
                   help='Replicates per model (default: 5).')
    p.add_argument('--models', default=','.join(MODELS),
                   help=f'Comma-separated model labels to run (default: all). '
                        f'Choices: {", ".join(MODELS)}.')
    p.add_argument('--results-root',
                   default=os.path.join(_ROOT, 'results', 'batches', 'hl_batches'),
                   help='Where to write <shot>_shot/ (default: results/batches/hl_batches).')
    p.add_argument('--print-prompts', action='store_true',
                   help='Print the system/user prompts for --shot and exit (no API calls).')
    args = p.parse_args(argv)

    labels = [m.strip() for m in args.models.split(',') if m.strip()]
    bad = [m for m in labels if m not in MODELS]
    if bad:
        raise SystemExit(f"Unknown model label(s) {bad}; choices: {list(MODELS)}")

    if args.print_prompts:
        if args.shot == 'zero':
            system, user = ZERO_SHOT_SYSTEM, ZERO_SHOT_USER
        elif args.shot == 'frag':
            system, user = FRAG_SHOT_SYSTEM, FRAG_SHOT_USER
        else:
            system, user = FEW_SHOT_SYSTEM, _few_shot_user()
        print(f'===== {args.shot}-shot SYSTEM =====\n{system}')
        print(f'===== {args.shot}-shot USER =====\n{user}')
        return 0

    run_batch(args.shot, args.replicates, labels, args.results_root)
    return 0


if __name__ == '__main__':
    sys.exit(main())
