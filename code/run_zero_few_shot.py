#!/usr/bin/env python3
"""
run_zero_few_shot.py - Non-agentic baseline: single-shot molecule generation
(no tools, no adversary, no multi-turn refinement) across all 5 proposers used
elsewhere in this study, for zero-shot, few-shot, and fragment-suggested
zero-shot prompting.

Ported from GPT_ANT_zero-shot.py / GPT_ANT_ONE_SHOT.py (prototypes that covered
only OpenAI + Anthropic, and whose few-shot prompt was wired to an unrelated
HOMO-LUMO-gap task) -- extended here to all 5 proposers (OpenAI, Anthropic,
Gemini, kimi-k2.6, deepseek-v4-pro, the last two via Ollama) and to the same
HMGCR docking task used throughout the rest of the study:

  zero-shot: system prompt only states the task; the user message names just
             the protein (no example molecules) -- verbatim from
             GPT_ANT_zero-shot.py's dock_task_specific_prompt / first_prompt.
  few-shot:  system prompt is docking_module.task_specific_prompt (the same
             text GPT_ANT_ONE_SHOT.py had defined locally but never wired up)
             plus that file's own "learn trends, then propose" scaffolding;
             the user message is the full code/adversarial_set.md molecule/
             docking-score list -- verbatim from GPT_ANT_ONE_SHOT.py.
  frag-shot: same task framing and bare-protein-name user message as
             zero-shot, but the system prompt also suggests a fixed set of
             base rings and functional groups to build molecules from.

Each (shot-type, model, replicate) is a single API call: no tools, no critique,
no loop. Output is written to match the existing per-batch convention exactly
(results/batches/<zero_shot|few_shot|frag_shot>/<model>/rep<N>/<...>.md + .json
sidecar + a manifest.json at the zero_shot/few_shot/frag_shot level) so
`analyze_replicates.py --batch-dir results/batches/zero_shot` (or few_shot /
frag_shot) works completely unmodified -- same real-Vina-redock analysis as
every agentic run in this study.

Usage:
  python3 code/run_zero_few_shot.py --shot zero --replicates 5
  python3 code/run_zero_few_shot.py --shot few --replicates 5
  python3 code/run_zero_few_shot.py --shot frag --replicates 5
  python3 code/run_zero_few_shot.py --shot zero --replicates 1 --models openai,gemini
"""
import os
import sys
import re
import json
import time
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root

from docking_module import task_specific_prompt  # noqa: E402  (the docking few-shot prompt)

_DEFAULT_API_TIMEOUT = 120.0

OPENAI_DEFAULT_MODEL = 'gpt-5.2'
ANTHROPIC_DEFAULT_MODEL = 'claude-haiku-4-5-20251001'
GEMINI_DEFAULT_MODEL = 'gemini-3-flash-preview'
OLLAMA_HOST = 'https://ollama.com'

# label -> (provider, default_model). 'ollama' proposers each get their own
# label (kimi-k2.6 / deepseek-v4-pro) since a single manifest here covers both,
# unlike the agentic batches which run one Ollama model per manifest.
MODELS = {
    'openai': ('openai', OPENAI_DEFAULT_MODEL),
    'anthropic': ('anthropic', ANTHROPIC_DEFAULT_MODEL),
    'gemini': ('gemini', GEMINI_DEFAULT_MODEL),
    'kimi-k2.6': ('ollama', 'kimi-k2.6'),
    'deepseek-v4-pro': ('ollama', 'deepseek-v4-pro'),
    # Additional Ollama-cloud-only proposers (no agentic-loop role in this study,
    # baselines only). API names verified to resolve via a direct ollama.chat()
    # smoke call before adding here; 'gpt-oss-*' labels differ from their model
    # tags (which contain ':') purely for a cleaner --models CLI value / filename.
    'gemma4': ('ollama', 'gemma4'),
    'glm-5.2': ('ollama', 'glm-5.2'),
    'nemotron-3-ultra': ('ollama', 'nemotron-3-ultra'),
    'nemotron-3-super': ('ollama', 'nemotron-3-super'),
    'nemotron-3-nano': ('ollama', 'nemotron-3-nano:30b'),
    'gpt-oss-20b': ('ollama', 'gpt-oss:20b'),
    'gpt-oss-120b': ('ollama', 'gpt-oss:120b'),
}


# --- Prompts -----------------------------------------------------------------
# Verbatim from GPT_ANT_zero-shot.py.
ZERO_SHOT_SYSTEM = '''# You are a drug design assistant. Your task is to design a new molecules
with the best possible docking score (the most negative) to a particular protein target, given in the first user message.
You will deliver up to five potential molecules in SMILES format, along with reasoning for why you chose those molecules
and an estimate of their docking scores.
'''
ZERO_SHOT_USER = 'HMGCR'

# zero-shot, but with a fixed set of suggested rings/functional groups to build
# from -- same task framing as ZERO_SHOT_SYSTEM, same bare-protein-name user
# message, just with fragment suggestions appended to the system prompt.
FRAG_SHOT_SYSTEM = '''# You are a drug design assistant. Your task is to design a new molecules
with the best possible docking score (the most negative) to a particular protein target, given in the first user message.
You will deliver up to five potential molecules in SMILES format, along with reasoning for why you chose those molecules
and an estimate of their docking scores.

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
- I
- C#N
- C(=O)O(C(C)C)
- C#C(SC)
- C(C(=O)[O-])
- C(C)
- C=C([N+](=O)[O-])
- C(N)
- C([O-])
- CC(N(C)C)
'''
FRAG_SHOT_USER = 'HMGCR'

# Verbatim structure from GPT_ANT_ONE_SHOT.py, with the docking prompt
# (task_specific_prompt) swapped in for the file's unused HL_task_specific_prompt.
FEW_SHOT_SYSTEM = f'''
{task_specific_prompt}

## You will first:
- Read the list of molecule SMILES and scores
- Ascertain any features of the molecules that contribute to the desired score. For example, if,
from one molecule to the next, the addition of an O group makes the score better.
- Gather all of these trends across all of the molecules.

## Once you have ascertained the trends:
- Use the trends you learned to suggest 1-5 new molecules that obey the trends you found
and which should have a better score than the molecules in the list.
- Provide reasoning as to why you created those new molecules.
- Estimate the new scores.
'''


def _few_shot_user() -> str:
    with open(os.path.join(_HERE, 'adversarial_set.md'), 'r') as f:
        context = f.read()
    return f'\n  Here is a list of molecules and their docking scores:\n  {context}\n'


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
# Same client-construction pattern as OpenAIActor/AnthropicActor/GeminiActor/
# OllamaAdversary (lazy import, max_retries=0 + timeout so a permanent 401/429
# fails fast instead of hanging on SDK backoff), but parameterized by an
# arbitrary system prompt rather than hardcoded to ADVERSARY_INSTRUCTIONS.

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


# Labels where think=True is impractical for a single-shot call: observed
# 2min-23min+ variable latency for kimi-k2.6 with no way to predict which,
# unacceptable for what's supposed to be a fast non-agentic baseline. Same
# fix already proven for kimi-k2.6 on the adversary/critique path (see
# OllamaAdversary.critique() in molopt.py). deepseek-v4-pro has been
# reliable with think=True (5/5 zero+few shot, no timeouts) so stays as-is.
# The 7 newer Ollama-only labels default to NO_THINK too: their think=True
# single-call latency hasn't been characterized, and a 3-shot x 5-rep batch
# (105 calls) is too large to risk an unpredictable multi-minute hang on an
# unverified model. Revisit per-label if think=True is specifically wanted.
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


# --- Output: .md + JSON sidecar matching the existing parse format ----------
# analyze_replicates.py (via verify_results.py's all_model_response_blocks)
# looks for a '# Initial model response:' section; the header line format
# ('# protein: ... | main model: ... | adversary: ...') is parsed for
# metadata but not required for SMILES extraction itself.

def write_run(results_dir, label, model, system, user, response, shot):
    os.makedirs(results_dir, exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    safe_model = model.replace(':', '').replace('/', '-')
    md_path = os.path.join(results_dir, f"{safe_model}_HMGCR_{timestamp}.md")
    sidecar_path = os.path.splitext(md_path)[0] + '.json'

    with open(md_path, 'w') as f:
        f.write(f'# {shot.capitalize()}-Shot Design Session - {timestamp}\n')
        f.write(f'# protein: HMGCR | main model: {model} (think=n/a) | adversary: none\n\n')
        f.write('# Initial model response:\n')
        f.write(response + '\n\n')
        f.write('# Session end: Done\n')

    messages = [
        {'role': 'system', 'content': system},
        {'role': 'user', 'content': user},
        {'role': 'assistant', 'content': response},
    ]
    payload = {
        'protein': 'HMGCR', 'shot': shot, 'label': label, 'model': model,
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
                'script': 'run_zero_few_shot.py', 'protein': 'HMGCR',
                'status': status, 'md_path': os.path.abspath(md_path),
                'sidecar_path': os.path.abspath(sidecar_path), 'exit_code': exit_code,
            })
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)

    manifest['batch_id'] = f'{shot}_shot'
    manifest['batch_dir'] = os.path.abspath(batch_dir)
    manifest['protein'] = 'HMGCR'
    manifest['replicates'] = replicates
    manifest['sets'] = labels
    manifest['models'] = {label: MODELS[label][1] for label in labels}
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\n{shot}-shot batch done. Manifest: {manifest_path}")
    print(f"Analyze with: fao-env/bin/python code/analyze_replicates.py --batch-dir {batch_dir}")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog='run_zero_few_shot.py',
        description='Non-agentic zero-shot / few-shot baseline across all 5 proposers.')
    p.add_argument('--shot', choices=['zero', 'few', 'frag'], required=True,
                   help='Which prompting mode to run.')
    p.add_argument('--replicates', type=int, default=5,
                   help='Replicates per model (default: 5).')
    p.add_argument('--models', default=','.join(MODELS),
                   help=f'Comma-separated model labels to run (default: all). '
                        f'Choices: {", ".join(MODELS)}.')
    p.add_argument('--results-root', default=os.path.join(_ROOT, 'results', 'batches'),
                   help='Where to write <shot>_shot/ (default: results/batches).')
    args = p.parse_args(argv)

    labels = [m.strip() for m in args.models.split(',') if m.strip()]
    bad = [m for m in labels if m not in MODELS]
    if bad:
        raise SystemExit(f"Unknown model label(s) {bad}; choices: {list(MODELS)}")

    run_batch(args.shot, args.replicates, labels, args.results_root)
    return 0


if __name__ == '__main__':
    sys.exit(main())
