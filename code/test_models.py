#!/usr/bin/env python3
"""
test_models.py - Smoke-test the candidate Ollama cloud models for molopt.py.

For each model listed in --models-file (default new_models_to_add.txt), checks:
  1. NAME     - the API model name resolves. We hypothesise dropping the
                '-cloud'/':cloud' suffix (e.g. 'gemma4:31b-cloud' -> 'gemma4:31b',
                'glm-5.2:cloud' -> 'glm-5.2'); the script tries that first, then
                the raw name as a fallback, and records which one worked.
  2. TOOLS    - the model emits a tool call when asked (mirrors molopt.py's
                ollama.chat(..., tools=...) tool-calling path).
  3. THINK    - the model can do think=True *and* tools in the same call. The
                old deepseek models could not; failures here flag a model for
                molopt.py's NO_THINK_MODELS set.
  4. ADVERSARY- one short turn: the main model proposes, the adversary critiques
                (uses the same OpenAI/Anthropic clients and instructions as
                molopt.py).

Only the ollama client is required. The adversary step uses the OpenAI/Anthropic
SDKs lazily, constructed exactly as in molopt.py, and keys are resolved the same
way (CLI flags, then env vars, then an optional .env file).
"""

import os
import sys
import time
import argparse

from ollama import Client as OllamaClient

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)  # repo root -- new_models_to_add.txt lives there
DEFAULT_MODELS_FILE = os.path.join(_ROOT, 'new_models_to_add.txt')


# --- Minimal .env loader (verbatim logic from molopt.py) --------------------

_ENV_KEYS = ('OLLAMA_API_KEY', 'OPENAI_API_KEY', 'ANTHROPIC_API_KEY')


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


# --- Name handling ----------------------------------------------------------

def strip_cloud(full: str) -> str:
    """Drop the trailing '-cloud' or ':cloud' suffix, then any dangling ':'.

    'gemma4:31b-cloud'    -> 'gemma4:31b'
    'qwen3.5:397b-cloud'  -> 'qwen3.5:397b'
    'glm-5.2:cloud'       -> 'glm-5.2'
    'kimi-k2.7-code:cloud'-> 'kimi-k2.7-code'
    """
    s = full
    for suf in ('-cloud', ':cloud'):
        if s.endswith(suf):
            s = s[: -len(suf)]
            break
    if s.endswith(':'):
        s = s[:-1]
    return s


def name_candidates(full: str):
    """Ordered API names to try: the stripped hypothesis first, then the raw name."""
    stripped = strip_cloud(full)
    cands = []
    if stripped and stripped != full:
        cands.append(stripped)
    cands.append(full)
    seen, out = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


# --- Trivial tool used to probe tool-calling --------------------------------
# A real chemistry tool (grow_cycle etc.) would invoke docking and make this
# smoke test slow. A trivial tool confirms the model emits tool_calls and that
# think+tools coexist - the actual tools are exercised by `molopt.py --self-test`.

def get_docking_score(smiles: str) -> str:
    """Return a placeholder docking score (kcal/mol) for a molecule given its
    SMILES string. Used here only to confirm the model will emit a tool call."""
    return f"docking score for {smiles}: -8.5 kcal/mol"


_TOOL = get_docking_score
_PROBE_SMILES = 'c1ccccc1'


def _is_model_not_found(err) -> bool:
    msg = str(err).lower()
    return any(k in msg for k in ('model not found', 'not found', '404', 'no such'))


# --- Adversary (verbatim from molopt.py) ------------------------------------

ADVERSARY_INSTRUCTIONS = '''
    You are a drug design assistant. You will recieve a proposal from  another model
    of novel molecules it has designed to bind to a particular protein target. The proposal will
    include reasoning as to why the model thinks those molecules will bind well, and estimated
    docking scores for each molecule. Your task is to analyze the proposal and find any flaws
    in the reasoning or estimation of the docking scores. You should then suggest modifications
    to the proposed molecules that would make them more likely to bind well, and provide reasoning
    for why those modifications would help.

    The other model has access to the following tools, and you may suggest that it use these tools to
    gather more information or test out modifications to the proposed molecules:

    - grow_cycle: starts with a molecule SMILES and adds substituents to it, docks them, and returns
                  a list of molecules and scores.

    - replace_groups: starts with a molecule SMILES and replaces specific groups in it with new groups,
                      returning a list of new molecules and scores.

    - make_random_list: this tool generates a list of substituents of specified length (num_items).

    - related: this tool generates a list of molecules that are structurally related to a given molecule,
                and may be useful for exploring the chemical space around promising molecules.

    - lipinski: this tool evaluates a list of molecules for their drug-likeness based on Lipinski's rule of five.

    - dock_and_get_interacting_residues: this tool docks a single molecule and returns the residues in the protein
                that it is interacting with, as well as the types of interactions.

    - calculate_SAS_and_NP: this tool calculates the synthetic accessibility score and the natural product likeness
                score; the former rates ease of synthesis on a 1 to 10 scale, and the latter rates similarity to natural
                products on a -5 to 5 scale.
'''


class OpenAIAdversary:
    """OpenAI Responses-API adversary (matches molopt.py)."""

    def __init__(self, model: str, api_key: str):
        from openai import OpenAI  # lazy import
        self.model = model
        self.client = OpenAI(api_key=api_key)

    def critique(self, prompt: str) -> str:
        resp = self.client.responses.create(
            model=self.model,
            instructions=ADVERSARY_INSTRUCTIONS,
            input=prompt,
        )
        return resp.output_text


class AnthropicAdversary:
    """Anthropic Messages-API adversary (matches molopt.py)."""

    def __init__(self, model: str, api_key: str, base_url: str = "https://api.anthropic.com"):
        from anthropic import Anthropic  # lazy import
        self.model = model
        # Explicit base_url: bypass any ambient ANTHROPIC_BASE_URL (matches molopt.py).
        self.client = Anthropic(api_key=api_key, base_url=base_url)

    def critique(self, prompt: str) -> str:
        resp = self.client.messages.create(
            model=self.model,
            system=ADVERSARY_INSTRUCTIONS,
            max_tokens=4096,
            messages=[{"role": "user", "content": prompt}],
        )
        parts = [b.text for b in resp.content if getattr(b, "type", None) == "text"]
        return "".join(parts).strip()


def make_adversary(provider: str, model: str, openai_key: str, anthropic_key: str,
                   anthropic_base_url: str = "https://api.anthropic.com"):
    if provider == 'openai':
        if not openai_key:
            return None
        return OpenAIAdversary(model, openai_key)
    if provider == 'anthropic':
        if not anthropic_key:
            return None
        return AnthropicAdversary(model, anthropic_key, base_url=anthropic_base_url)
    return None


# --- Individual checks ------------------------------------------------------

def check_name(ollama, full):
    """Try candidate API names; return (working_name, error_or_None)."""
    last_err = None
    for name in name_candidates(full):
        try:
            ollama.chat(
                model=name,
                messages=[{'role': 'user', 'content': 'Reply with the single word: OK'}],
                think=False,
                options={'num_predict': 16},
            )
            return name, None
        except Exception as err:
            last_err = err
            if _is_model_not_found(err):
                continue  # try the next candidate name
            return None, err  # other error (auth/network) - stop here
    return None, last_err


def check_tools(ollama, model, *, think):
    """Ask the model to call get_docking_score (mirrors molopt.py's tool path).

    A single round-trip: the tool is never executed and re-sent, so this cannot
    loop. num_predict caps the (possibly long) thinking trace so a verbose
    model can't drag one call out.
    """
    try:
        response = ollama.chat(
            model=model,
            messages=[{
                'role': 'user',
                'content': (f'What is the docking score for the SMILES {_PROBE_SMILES}? '
                            f'You MUST call the get_docking_score tool to answer.'),
            }],
            tools=[_TOOL],
            think=think,
            options={'num_predict': 512},
        )
        tool_calls = getattr(response.message, 'tool_calls', None) or []
        if tool_calls:
            names = [getattr(tc.function, 'name', '?') for tc in tool_calls]
            thinking = getattr(response.message, 'thinking', None)
            th = ' +thinking' if (think and thinking) else ''
            return True, f'called {names}{th}'
        content = (getattr(response.message, 'content', '') or '').strip()
        snippet = (content[:80] + '…') if len(content) > 80 else content
        return False, f'no tool call; said: {snippet!r}'
    except Exception as err:
        return False, f'error: {err}'


def check_adversary(ollama, model, adversary):
    """One short turn: main model proposes two molecules, adversary critiques."""
    if adversary is None:
        return 'skipped (no key / --no-adversary)'
    try:
        prop = ollama.chat(
            model=model,
            messages=[{'role': 'user', 'content':
                'Propose two simple drug-like molecules as SMILES with a one-line '
                'rationale each. Keep it under 80 words.'}],
            tools=[], think=False, options={'num_predict': 300},
        )
        proposal = (getattr(prop.message, 'content', '') or '').strip()
        if not proposal:
            return 'skipped (main model returned no proposal text)'
        crit = adversary.critique(proposal)
        if not crit:
            return 'adversary replied empty'
        return f'ok ({len(crit)} chars): {crit[:120]!r}'
    except Exception as err:
        return f'error: {err}'


# --- Driver -----------------------------------------------------------------

def load_models(path):
    models = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                models.append(line)
    return models


def build_arg_parser():
    p = argparse.ArgumentParser(description='Smoke-test candidate Ollama cloud models.')
    p.add_argument('--models-file', default=DEFAULT_MODELS_FILE)
    p.add_argument('--only', default=None, help='Comma-separated subset of model names to test.')
    # Ollama (mirrors molopt.py defaults).
    p.add_argument('--ollama-host', default='https://ollama.com', help='Ollama host (default: https://ollama.com).')
    p.add_argument('--ollama-key', default=None, help='Ollama bearer token (or env OLLAMA_API_KEY).')
    p.add_argument('--timeout', type=float, default=180.0,
                   help='Per-call HTTP timeout in seconds (default: 180). Bounds a '
                        'slow cloud cold-start so one stuck call cannot hang the run.')
    # Adversary (mirrors molopt.py defaults: openai / gpt-5.2).
    p.add_argument('--no-adversary', action='store_true')
    p.add_argument('--adversary', choices=['openai', 'anthropic'], default='openai',
                   help='Adversary provider (default: openai).')
    p.add_argument('--adversary-model', default=None,
                   help='Adversary model (default: gpt-5.2 for openai, claude-haiku-4-5-20251001 for anthropic).')
    p.add_argument('--openai-key', default=None, help='OpenAI API key (or env OPENAI_API_KEY).')
    p.add_argument('--anthropic-key', default=None, help='Anthropic API key (or env ANTHROPIC_API_KEY).')
    p.add_argument('--anthropic-base-url', default=None,
                   help='Anthropic API base URL (default: https://api.anthropic.com; bypasses ANTHROPIC_BASE_URL env).')
    return p


def main(argv=None):
    load_dotenv()
    args = build_arg_parser().parse_args(argv)

    # Ollama client - constructed exactly as in molopt.py.
    ollama_key = (args.ollama_key or os.environ.get('OLLAMA_API_KEY')
                  or os.environ.get('OLLAMA_KEY') or '')
    headers = {'Authorization': f'Bearer {ollama_key}'} if ollama_key else {}
    ollama = OllamaClient(host=args.ollama_host, headers=headers, timeout=args.timeout)

    # Adversary keys - resolved exactly as in molopt.py.
    openai_key = args.openai_key or os.environ.get('OPENAI_API_KEY') or ''
    anthropic_key = args.anthropic_key or os.environ.get('ANTHROPIC_API_KEY') or ''

    adv_model = args.adversary_model
    if adv_model is None:
        adv_model = 'gpt-5.2' if args.adversary == 'openai' else 'claude-haiku-4-5-20251001'
    adversary = None if args.no_adversary else make_adversary(
        args.adversary, adv_model, openai_key, anthropic_key,
        anthropic_base_url=args.anthropic_base_url or "https://api.anthropic.com")

    models = load_models(args.models_file)
    if args.only:
        wanted = {s.strip() for s in args.only.split(',') if s.strip()}
        models = [m for m in models if m in wanted or strip_cloud(m) in wanted]

    print(f'Ollama host : {args.ollama_host}  (key: {"yes" if ollama_key else "NO"})')
    print(f'Adversary   : {("skipped" if adversary is None else args.adversary + "/" + str(adv_model))}')
    print(f'Models      : {", ".join(models)}')
    print('=' * 78)

    results = []
    no_think_recs = []
    for full in models:
        t0 = time.time()
        print(f'\n## {full}')
        print(f'   candidates: {name_candidates(full)}')

        name, err = check_name(ollama, full)
        if name is None:
            print(f'   NAME      : FAIL - could not resolve ({err})')
            results.append((full, None, False, False, 'name-failed'))
            continue
        print(f'   NAME      : OK  -> {name}  (api name)')

        ok_tools, det = check_tools(ollama, name, think=False)
        print(f'   TOOLS     : {"OK  " if ok_tools else "FAIL"} - {det}')

        ok_think, det2 = check_tools(ollama, name, think=True)
        print(f'   THINK+TOOL: {"OK  " if ok_think else "FAIL"} - {det2}')
        if not ok_think:
            no_think_recs.append(name)

        adv = check_adversary(ollama, name, adversary)
        print(f'   ADVERSARY : {adv}')

        print(f'   ({time.time() - t0:.1f}s)')
        results.append((full, name, ok_tools, ok_think, adv))

    print('\n' + '=' * 78)
    print('SUMMARY')
    print(f'{"model":28} {"api name":22} tools think adv')
    for full, name, ok_t, ok_th, adv in results:
        if name is None:
            print(f'{full:28} {"-":22} -     -     name-failed')
            continue
        t = 'ok' if ok_t else 'NO'
        th = 'ok' if ok_th else 'NO'
        a = 'ok' if adv.startswith('ok') else ('skip' if adv.startswith('skipped') else 'NO')
        print(f'{full:28} {name:22} {t:5} {th:5} {a}')

    if no_think_recs:
        print('\nRecommend adding to molopt.py NO_THINK_MODELS (think+tools failed):')
        print('  ' + ', '.join(repr(m) for m in no_think_recs))
    else:
        print('\nAll models supported think=True + tools together.')


if __name__ == '__main__':
    main()