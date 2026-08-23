#!/usr/bin/env python3
"""
molopt.py - Command-line adversarial molecule optimization.

Refactor of Ollama_MolOpt.ipynb. Runs the same headless loop locally:
  1. An Ollama main model reasons over a molecule/docking-score list and may
     call chemistry tools (grow_cycle, replace_groups, make_random_list,
     related, lipinski, dock_and_get_interacting_residues, calculate_SAS_and_NP).
  2. An adversary model (OpenAI or Anthropic) critiques each proposal.
  3. The main model refines until it replies "Done" (or --max-turns is hit).

Every step is appended to a timestamped Markdown file under --results-dir.

Secrets are read from CLI flags or environment variables (optionally a .env
file in the working directory). No Colab / gradio dependencies.
"""

import os
import sys
import re
import json
import time
import argparse

# --- Path + NumPy-2 shim setup (must run before importing the helpers) ------
# The helpers live in ./code; put that dir on the path.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, 'code'))

# ODDT 0.7 calls np.in1d, removed in NumPy 2.x. np.isin is a drop-in.
# docking_module.py also shims this, but we set it first to be safe regardless
# of import order.
import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin

# --- Helper imports (tools, shared scoring_args, prompt templates) ----------
from MolPropOp import (  # noqa: F401  (imported names are used as tools / below)
    grow_cycle,
    replace_groups,
    make_random_list,
    related,
    lipinski,
)
from docking_module import (  # noqa: F401
    dock_and_get_interacting_residues,
    calculate_SAS_and_NP,
    scoring_args,
    task_specific_prompt,
    task_specific_tools,
)
import mock_tools

# Silence RDKit's "SMILES Parse Error" stderr noise. Invalid substituents the
# model passes are now reported back to it as 'invalid SMILES, skipped' entries
# in the tool results (see MolPropOp.py), so the stderr warnings are redundant.
# Re-enable them with --rdkit-verbose.
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Ollama client is always needed; import up front.
from ollama import Client as OllamaClient

# Models offered via the Ollama cloud endpoint, kept for --list-models
# convenience. These are the API names: the cloud catalog lists each with a
# '-cloud'/':cloud' suffix (e.g. 'gemma4:31b-cloud', 'glm-5.2:cloud'), but that
# suffix must be dropped for the chat API (verified by test_models.py).
OLLAMA_MODELS = [
    'gemma4:31b', 'glm-5.2', 'kimi-k2.7-code',
    'deepseek-v4-pro', 'qwen3.5:397b',
]

# Models for which Ollama "think" mode is auto-disabled. The previous deepseek
# generations could not do think + tool-calling simultaneously; deepseek-v4-pro
# can (verified by test_models.py), so this set is now empty. Add a model here
# only if think=True breaks its tool-calling.
NO_THINK_MODELS = set()

# --- Previous model set (kept for backwards compatibility) -----------------
# The notebook's original Ollama model list and its think-disabled subset.
# Uncomment to restore the old behavior (e.g. to re-run an older model).
#
# OLLAMA_MODELS = [
#     'deepseek-v3.1:671b', 'gpt-oss:120b', 'gpt-oss:20b',
#     'devstral-2:123b', 'cogito-2.1:671b',
#     'nemotron-3-nano:30b', 'gemini-3-flash-preview',
#     'kimi-k2:1t', 'kimi-k2.5', 'deepseek-v3.2',
# ]
# NO_THINK_MODELS = {'deepseek-v3.1:671b', 'deepseek-v3.2'}

# Callables the main model may invoke, in stable order.
TOOL_FUNCTIONS = [
    grow_cycle,
    replace_groups,
    make_random_list,
    related,
    lipinski,
    dock_and_get_interacting_residues,
    calculate_SAS_and_NP,
]
AVAILABLE_FUNCTIONS = {fn.__name__: fn for fn in TOOL_FUNCTIONS}

# --- Prompt text (ported verbatim from notebook cell 11) --------------------

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


def build_system_message() -> str:
    """Build the main model's system message (notebook cell 11 `sys_message`).

    Uses the imported task_specific_prompt / task_specific_tools templates.
    """
    return f'''
{task_specific_prompt}

## You will first:
- Read the list of molecule SMILES and scores
- Ascertain any features of the molecules that contribute to the desired score. For example, if,
from one molecule to the next, the addition of an O group makes the score better.
- Gather all of these trends across all of the molecules.

## If you need additional information to ascertain the trends, such as more modified
molecules and their docking scores, you have tools you can call to generate new
molecules and get their docking scores. You can use these tools as many times as you want
to gather information on the trends. *NOTE: if you choose to add a phenyl group to a molecule,
use the SMILES 'c7ccccc7', so that it does not interfere with other rings in the molecule that
may already use numbers 1-6 in their SMILES notation.

The tools you have available include:

- grow_cycle: starts with a molecule SMILES and adds substituents to it, docks them, and returns
              a list of molecules and scores. You can use this tool to further explore modifications
              to promising molecules that you find in the input data. You can provide a list of
              substituents to add, or use the predefined sets: e_withdraw (electron withdrawing),
              e_donate (electron donating), withdraw_with_linkers (electron withdrawing with linkers),
              donate_with_linkers (electron donating with linkers). You can also generate a random list
              of substituents with the make_random_list tool and use that as input to grow_cycle.

- replace_groups: starts with a molecule SMILES and replaces specific groups in it with new groups, returning a list of new
                  molecules and scores. This tool allows you to test specific hypotheses about how replacing certain
                  groups in a molecule might affect binding affinity. You can specify which groups to replace and
                  what to replace them with, or use the predefined sets of substituents mentioned above. You can also
                  generate a random list of substituents with the make_random_list tool and use that as input to replace_groups.

- make_random_list: this tool generates a list of substituents of specified length (num_items). It draws from the predefined lists:
                    e_withdraw (electron withdrawing), e_donate (electron donating), withdraw_with_linkers
                    (electron withdrawing with linkers), donate_with_linkers (electron donating with linkers).
                    Use this tool when you want to get a broad sense of how different modifications affect binding affinity,
                    without having a specific hypothesis in mind.

- related: this tool generates a list of molecules that are structurally related to a given molecule, and
           may be useful for exploring the chemical space around promising molecules you find in the input
           data. It returns a list of related molecules and a few properties.

- lipinski: this tool evaluates a list of molecules for their drug-likeness based on Lipinski's rule of five,
            which is a set of guidelines for determining whether a molecule is likely to be an orally active
            drug in humans. This tool can help you ensure that the molecules you are proposing not only have
            good docking scores but also have properties that make them more likely to be successful as drugs.
            QED (quantitative estimate of drug-likeness) is a score between 0 and 1 that summarizes how
            drug-like a molecule is, with 1 being the most drug-like. A higher QED score indicates that a
            molecule has properties that are more consistent with known drugs, such as appropriate molecular
            weight, lipophilicity, and number of hydrogen bond donors and acceptors.

{task_specific_tools}

## Once you have ascertained the trends:
- Use the trends you learned to suggest 1-5 new molecules that obey the trends you found
and which should have a better score than the molecules in the list.
- Provide reasoning as to why you created those new molecules.
- Estimate the new scores.

## You may ask the user for clarification if needed, but try to use the tools to gather as much information as you
can before asking for clarification.

## In further turns, you will also receive feedback from an adversary model that is trying to find flaws
in your reasoning and suggest improvements to your proposed molecules. You should use this feedback to
refine your understanding of the trends, run new experiments with the tools to gather more information,
and improve your proposed molecules in subsequent turns.

## If you have identified good potential hits, evaluate the Lipinski properties of the proposed molecules
and use that information to further refine your proposals, keeping in mind that you want to propose molecules
that not only have good docking scores but also have good drug-like properties.

## Once you have reached a point where you think you have proposed the best possible molecules based on the
trends, tool results and the adversary feedback, reply with only one word: "Done". This will signal that you have
finished the task and will not propose any more molecules.
'''


# --- Adversary abstraction --------------------------------------------------

# --- API-error classification ------------------------------------------------
#
# Both OpenAI and Anthropic SDKs default to max_retries=2 with exponential
# backoff and a 600s per-call timeout. For a PERMANENT error -- bad auth (401/
# 403) or an out-of-credits 429 -- that backoff makes the call appear to hang
# for many minutes at 0% CPU instead of failing fast. So we build both
# adversary clients with max_retries=0 + an explicit timeout (chat_turn's
# manual loop handles retries for the Ollama proposer), and refuse to retry
# errors that can never succeed. Mirrors molopt_oa.py's fail-fast pattern.

_DEFAULT_API_TIMEOUT = 120.0  # seconds per call; well under the SDK's 600s default


def _is_retryable(err) -> bool:
    """True for transient errors worth retrying; False for permanent ones.

    Permanent (do not retry): 401/403 (auth), and a 429 that is actually
    billing/quota exhaustion (e.g. 'credit_balance_exhausted' /
    'insufficient_quota') rather than a transient rate limit. Connection
    errors and 5xx have no status_code and are treated as retryable.
    """
    sc = getattr(err, 'status_code', None)
    if sc in (401, 403):
        return False
    if sc == 429:
        low = str(err).lower()
        if any(k in low for k in ('credit', 'quota', 'billing', 'insufficient')):
            return False
    return True


class OpenAIAdversary:
    """OpenAI Responses-API adversary (matches notebook cell 11)."""

    def __init__(self, model: str, api_key: str, timeout: float = _DEFAULT_API_TIMEOUT):
        from openai import OpenAI  # lazy import
        self.model = model
        # max_retries=0: avoid the SDK's minutes-long backoff on a permanent 429/401.
        self.client = OpenAI(api_key=api_key, timeout=timeout, max_retries=0)

    def critique(self, prompt: str) -> str:
        resp = self.client.responses.create(
            model=self.model,
            instructions=ADVERSARY_INSTRUCTIONS,
            input=prompt,
        )
        return resp.output_text


class AnthropicAdversary:
    """Anthropic Messages-API adversary."""

    def __init__(self, model: str, api_key: str, base_url: str = "https://api.anthropic.com",
                 timeout: float = _DEFAULT_API_TIMEOUT):
        from anthropic import Anthropic  # lazy import
        self.model = model
        # Explicit base_url: bypass any ambient ANTHROPIC_BASE_URL (e.g. a local
        # proxy that doesn't serve Anthropic model names), routing to the real API.
        # max_retries=0 + timeout: avoid the SDK's minutes-long backoff on a permanent 429/401.
        self.client = Anthropic(api_key=api_key, base_url=base_url,
                                timeout=timeout, max_retries=0)

    def critique(self, prompt: str) -> str:
        resp = self.client.messages.create(
            model=self.model,
            system=ADVERSARY_INSTRUCTIONS,
            max_tokens=4096,
            messages=[{"role": "user", "content": prompt}],
        )
        # Extract text from the content blocks.
        parts = [b.text for b in resp.content if getattr(b, "type", None) == "text"]
        return "".join(parts).strip()


def make_adversary(provider: str, model: str, openai_key: str, anthropic_key: str,
                   anthropic_base_url: str = "https://api.anthropic.com",
                   timeout: float = _DEFAULT_API_TIMEOUT):
    if provider == 'openai':
        if not openai_key:
            raise SystemExit("OpenAI adversary needs --openai-key or OPENAI_API_KEY.")
        return OpenAIAdversary(model, openai_key, timeout=timeout)
    if provider == 'anthropic':
        if not anthropic_key:
            raise SystemExit("Anthropic adversary needs --anthropic-key or ANTHROPIC_API_KEY.")
        return AnthropicAdversary(model, anthropic_key, base_url=anthropic_base_url,
                                  timeout=timeout)
    raise SystemExit(f"Unknown adversary provider: {provider!r} (use 'openai' or 'anthropic').")


# --- Minimal .env loader (no extra dependency) ------------------------------

_ENV_KEYS = ('OLLAMA_API_KEY', 'OPENAI_API_KEY', 'ANTHROPIC_API_KEY')


def load_dotenv(path: str = '.env') -> None:
    """Populate os.environ from a simple .env file if present.

    Only sets keys that are not already defined in the environment, so real
    env vars always win. Handles KEY=value and KEY="value" lines; ignores
    blanks and # comments.
    """
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


# --- Main-model tool-calling turn -------------------------------------------

def _vprint(verbose: bool, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)


def _msg_role(msg) -> str:
    if isinstance(msg, dict):
        return msg.get('role', '')
    return getattr(msg, 'role', '')


def _msg_content(msg) -> str:
    if isinstance(msg, dict):
        return msg.get('content', '') or ''
    return getattr(msg, 'content', '') or ''


def last_assistant_text(messages) -> str:
    """Most recent non-empty assistant message text.

    Ollama sometimes ends a turn with an empty-content assistant message (the
    real text was in an earlier message alongside a tool call). Returning the
    literal last message would hand the adversary empty input. This scans back
    for the latest assistant message that actually has text.
    """
    for msg in reversed(messages):
        if _msg_role(msg) == 'assistant':
            content = _msg_content(msg)
            if content and content.strip():
                return content
    return ''


# --- Messages sidecar (for --resume) ----------------------------------------

def _msg_to_serializable(msg) -> dict:
    """Convert a messages-list entry (dict or ollama Message) to a JSON-safe dict."""
    if isinstance(msg, dict):
        return msg
    # ollama's Message is a pydantic model; model_dump() yields plain JSON types.
    if hasattr(msg, 'model_dump'):
        return msg.model_dump()
    return {'role': getattr(msg, 'role', None), 'content': getattr(msg, 'content', None)}


def _serialize_messages(messages) -> list:
    return [_msg_to_serializable(m) for m in messages]


def _strip_for_api(messages) -> list:
    """Prepare loaded sidecar messages for sending back to ollama.chat.

    Drops the assistant `thinking`/`images` fields (model-internal; not expected
    in inputs) while keeping `tool_calls`/`tool_name` so the tool-calling
    conversation round-trips correctly.
    """
    out = []
    for m in messages:
        m = _msg_to_serializable(m)
        m = dict(m)
        m.pop('thinking', None)
        m.pop('images', None)
        out.append(m)
    return out


def write_sidecar(sidecar_path: str, *, args, messages, written_at_turn, status) -> None:
    """Atomically write the messages list + run metadata as a JSON sidecar.

    Written incrementally (after each turn) so a run killed mid-way still
    leaves a usable sidecar of the last completed turn. Atomic via tmp+rename
    so a kill mid-write cannot leave a half-written file.
    """
    payload = {
        'protein': args.protein,
        'model': args.model,
        'think': args.think,
        'adversary': args.adversary,
        'adversary_model': args.adversary_model,
        'max_turns': args.max_turns,
        'max_tool_calls': args.max_tool_calls,
        'written_at_turn': written_at_turn,
        'status': status,
        'messages': _serialize_messages(messages),
    }
    tmp = sidecar_path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(payload, f, indent=2, default=str)
    os.replace(tmp, sidecar_path)


def load_sidecar(path: str) -> list:
    """Load the messages list from a JSON sidecar, ready to feed to ollama.chat.

    Ensures the first message is the system prompt; if the sidecar is missing
    one (e.g. an older/edited file), prepend the current system message.
    """
    with open(path, 'r') as f:
        payload = json.load(f)
    messages = _strip_for_api(payload.get('messages', []))
    if not messages or _msg_role(messages[0]) != 'system':
        messages.insert(0, {'role': 'system', 'content': build_system_message()})
    return messages


_DONE_RE = re.compile(r'Done[.!?]*\s*$', re.IGNORECASE)


def is_done(text: str) -> bool:
    """True if the model signalled completion with "Done".

    The system prompt asks the model to reply with only the word "Done" when
    finished. In practice models often append "Done." (with punctuation) at the
    end of a final proposal, or wrap it on its own last line. Treat the response
    as Done if the whole text is just "Done" (±punctuation) or its final line
    is just "Done" (±punctuation). A bare "Done" mid-paragraph does not count.
    """
    if not text or not text.strip():
        return False
    s = text.strip()
    if _DONE_RE.fullmatch(s):
        return True
    last_line = s.splitlines()[-1].strip()
    return bool(_DONE_RE.fullmatch(last_line))


def chat_turn(messages, prompt, *, ollama, model, think, max_retries, max_tool_calls, verbose):
    """Send one user prompt and let the main model run until it stops calling tools.

    Returns (messages, last_assistant_content). Replaces the notebook's
    `chat_turn`: adds connection-error retry, tool-error isolation, and a
    per-turn tool-call cap so a stuck model (e.g. repeatedly retrying invalid
    inputs it can't see failed) can't loop forever.
    """
    messages = list(messages)
    messages.append({'role': 'user', 'content': prompt})

    trace = []  # captured thinking/content/tool-call lines, returned for --trace
    tool_rounds = 0
    while True:
        # Once the cap is hit, drop the tools so the model must emit text.
        force_tools = [] if tool_rounds >= max_tool_calls else TOOL_FUNCTIONS
        if tool_rounds >= max_tool_calls:
            _vprint(verbose, f"  [tool-call cap ({max_tool_calls}) reached; "
                             f"forcing a text-only response]")

        response = None
        last_err = None
        for attempt in range(1, max_retries + 1):
            try:
                response = ollama.chat(
                    model=model,
                    messages=messages,
                    tools=force_tools,
                    think=think,
                )
                break
            except Exception as err:  # connection / transient errors
                last_err = err
                _vprint(verbose, f"  [ollama error, attempt {attempt}/{max_retries}: {err}]")
                if not _is_retryable(err):
                    raise  # permanent (auth / out-of-credits): fail fast, don't loop
                if attempt < max_retries:
                    time.sleep(2 * attempt)
        if response is None:
            # Surface the error to the model so it can adapt; keep the loop alive.
            msg = (f"The previous call failed with a connection error: {last_err}. "
                   f"Please proceed from the last step.")
            _vprint(verbose, f"  [giving up after {max_retries} retries; asking model to continue]")
            messages.append({'role': 'user', 'content': msg})
            continue

        messages.append(response.message)
        _vprint(verbose, '-' * 72)
        _vprint(verbose, "Thinking: ", getattr(response.message, 'thinking', None))
        _vprint(verbose, '-' * 72)
        _vprint(verbose, "Content: ", response.message.content)
        _vprint(verbose, '-' * 72)
        trace.append(f"### thinking\n{getattr(response.message, 'thinking', None) or ''}")
        trace.append(f"### content\n{response.message.content or ''}")

        tool_calls = getattr(response.message, 'tool_calls', None)
        if tool_calls and tool_rounds < max_tool_calls:
            tool_rounds += 1
            for tc in tool_calls:
                fn_name = tc.function.name
                fn_args = tc.function.arguments
                _vprint(verbose, f"Calling {fn_name} with arguments {fn_args}")
                if fn_name in AVAILABLE_FUNCTIONS:
                    try:
                        result = AVAILABLE_FUNCTIONS[fn_name](**fn_args)
                    except Exception as err:
                        # Feed the error back so the model can self-correct
                        # (e.g. a typo'd argument name), instead of crashing.
                        result = f"Tool '{fn_name}' raised an error: {err}. " \
                                 f"Check the argument names/types and retry."
                        _vprint(verbose, f"  [tool error: {err}]")
                    _vprint(verbose, f"Result: {result}")
                    _vprint(verbose, '-' * 72)
                    trace.append(f"### tool call: {fn_name}\nargs: {fn_args}\nresult: {result}")
                    messages.append({'role': 'tool', 'tool_name': fn_name, 'content': str(result)})
                else:
                    _vprint(verbose, f"  [unknown tool {fn_name}; ignored]")
                    messages.append({'role': 'tool', 'tool_name': fn_name,
                                     'content': f"Error: unknown tool '{fn_name}'."})
        else:
            # Either the model stopped calling tools, or it kept calling them
            # past the cap. In the latter case the model may have produced no
            # usable text (its last message is just a tool request), so ask
            # once for a text-only summary with no tools available, then stop.
            if tool_calls:
                _vprint(verbose, "  [over cap; requesting a text summary without tools]")
                messages.append({'role': 'user', 'content':
                    "You have reached the tool-call limit for this turn. Do not call any "
                    "more tools. Summarize your best proposed molecules so far, their "
                    "estimated docking scores, and your reasoning. Do NOT say 'Done' — "
                    "you will receive adversary feedback next and then refine your proposals."})
                try:
                    final = ollama.chat(model=model, messages=messages, tools=[], think=think)
                    messages.append(final.message)
                    _vprint(verbose, '-' * 72)
                    _vprint(verbose, "Content: ", getattr(final.message, 'content', None))
                    _vprint(verbose, '-' * 72)
                    trace.append(f"### content (over-cap summary)\n{getattr(final.message, 'content', None) or ''}")
                except Exception as err:
                    _vprint(verbose, f"  [final summary call failed: {err}]")
            break  # turn is done

    last_content = last_assistant_text(messages)
    return messages, last_content, trace


# --- Session runner ---------------------------------------------------------

def run_session(args) -> str:
    """Run the full adversarial session. Returns the results file path."""
    # Configure the shared scoring state (same mutable object the helpers use).
    scoring_args[0] = os.cpu_count()
    scoring_args[1] = args.protein

    ollama_key = args.ollama_key or os.environ.get('OLLAMA_API_KEY') or os.environ.get('OLLAMA_KEY') or ''
    headers = {'Authorization': f'Bearer {ollama_key}'} if ollama_key else {}
    # timeout passed to the ollama Client -> httpx (bounds a hung connection).
    ollama = OllamaClient(host=args.ollama_host, headers=headers, timeout=args.api_timeout)

    adversary = make_adversary(
        args.adversary, args.adversary_model,
        args.openai_key or os.environ.get('OPENAI_API_KEY') or '',
        args.anthropic_key or os.environ.get('ANTHROPIC_API_KEY') or '',
        anthropic_base_url=args.anthropic_base_url or "https://api.anthropic.com",
        timeout=args.api_timeout,
    )

    think = args.think  # already resolved (None -> auto) in main()

    os.makedirs(args.results_dir, exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    safe_model = args.model.replace(':', '').replace('/', '-')
    results_path = os.path.join(args.results_dir, f"{safe_model}_{args.protein}_{timestamp}.md")
    # JSON messages sidecar lives next to the results .md (same stem). Always-on;
    # written incrementally so a stopped/killed run can be resumed via --resume.
    sidecar_path = os.path.splitext(results_path)[0] + '.json'
    with open(results_path, 'w') as f:
        f.write(f'# Adversarial Design Session - {timestamp}\n')
        f.write(f'# protein: {args.protein} | main model: {args.model} (think={think}) '
                f'| adversary: {args.adversary}/{args.adversary_model}\n')
        if args.resume:
            f.write(f'# resumed from: {args.resume}\n')
        f.write('\n')

    def log(section, text):
        with open(results_path, 'a') as f:
            f.write(f'\n{section}\n{text}\n')

    if args.resume:
        # Resume: seed the conversation from a prior run's JSON sidecar instead
        # of the context file, and skip the initial model turn. We pick up from
        # the last assistant message already in the sidecar.
        print(f"Resuming from sidecar -> {args.resume}")
        messages = load_sidecar(args.resume)
        last = last_assistant_text(messages)
        log('# Resumed from sidecar:', args.resume)
        log('# Last assistant text at resume:', last or '(none)')
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=0, status='resumed')
        turn = 0
    else:
        messages = [{'role': 'system', 'content': build_system_message()}]

        # Initial prompt: the molecule / docking-score list for this protein.
        with open(args.context_file, 'r') as f:
            context = f.read()
        first_prompt = f'\n  Here is a list of molecules and their docking scores:\n  {context}\n'

        print(f"Starting session -> {results_path}")
        messages, last, trace = chat_turn(
            messages, first_prompt,
            ollama=ollama, model=args.model, think=think,
            max_retries=args.max_retries, max_tool_calls=args.max_tool_calls,
            verbose=not args.quiet,
        )
        log('# Initial model response:', last)
        if args.trace:
            log('# Trace:', '\n'.join(trace))
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=0, status='in_progress')
        turn = 0

    while not is_done(last) and turn < args.max_turns:
        turn += 1
        print(f"\n=== Turn {turn}/{args.max_turns} ===")
        if last and last.strip():
            try:
                adv = adversary.critique(last)
                log('# Adversary feedback:', adv)
                _vprint(not args.quiet, f"[adversary {args.adversary}/{args.adversary_model} "
                                        f"replied ({len(adv)} chars)]")
            except Exception as err:
                # Don't let a single adversary API error kill a long run.
                adv = (f"The adversary model could not be reached (error: {err}). "
                       f"Review your latest proposal yourself, correct any flaws you "
                       f"can identify, and present your best molecules with estimated "
                       f"scores. Say 'Done' if you are finished.")
                log('# Adversary feedback: [unavailable]', adv)
                _vprint(not args.quiet, f"[adversary error: {err}]")
        else:
            # Model ended its turn with no text; nudge it to summarize instead
            # of calling the adversary with empty input (which OpenAI rejects).
            adv = ("Your last response had no text. Summarize your best proposed "
                   "molecules, their estimated docking scores, and your reasoning. "
                   "Say 'Done' if you are finished.")
            log('# Adversary feedback: [skipped - empty model response]', adv)
            _vprint(not args.quiet, "[adversary skipped - empty model response]")

        messages, last, turn_trace = chat_turn(
            messages, adv,
            ollama=ollama, model=args.model, think=think,
            max_retries=args.max_retries, max_tool_calls=args.max_tool_calls,
            verbose=not args.quiet,
        )
        log('# Model response:', last)
        if args.trace:
            log('# Trace:', '\n'.join(turn_trace))
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=turn, status='in_progress')

    status = 'Done' if is_done(last) else f'MAX_TURNS_REACHED (last={last!r})'
    with open(results_path, 'a') as f:
        f.write(f'\n# Session end: {status}\n')
    write_sidecar(sidecar_path, args=args, messages=messages,
                  written_at_turn=turn, status='Done' if is_done(last) else 'max_turns_reached')
    print(f"\nSession ended ({status}). Results: {results_path}")
    print(f"Messages sidecar: {sidecar_path}")
    return results_path


# --- Self-test (no LLM keys required) ---------------------------------------

def self_test() -> None:
    """Exercise the chemistry tools directly to confirm the stack works."""
    print("self-test: importing tools and running chemistry calls...")
    print("\n-- dock_and_get_interacting_residues('c1ccc(O)cc1') --")
    print(dock_and_get_interacting_residues('c1ccc(O)cc1'))
    print("\n-- calculate_SAS_and_NP(['c1ccccc1', 'c1ccc(O)cc1']) --")
    print(calculate_SAS_and_NP(['c1ccccc1', 'c1ccc(O)cc1']))
    print("\nself-test OK.")


# --- CLI --------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='molopt.py',
        description='Adversarial molecule optimization: Ollama main model + OpenAI/Anthropic adversary.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Show the built-in Ollama model set, then run with the default model
  python3 molopt.py --list-models
  python3 molopt.py --protein HMGCR

  # Specific cloud model (API name, no '-cloud' suffix) with think forced on
  python3 molopt.py --model gemma4:31b --think

  # Anthropic adversary instead of the default OpenAI one
  python3 molopt.py --model qwen3.5:397b --adversary anthropic --adversary-model claude-haiku-4-5-20251001

  # Quick chemistry-stack check with no LLM keys
  python3 molopt.py --self-test

Keys are read from the environment (OLLAMA_API_KEY, OPENAI_API_KEY,
ANTHROPIC_API_KEY) or a .env file in the working directory; CLI flags override.
Source ~/.zshrc first if the keys live there.
""",
    )
    p.add_argument('--protein', default='HMGCR',
                   help='Docking target, a dockstring target name (default: HMGCR). '
                        'Docking works for any of dockstring\'s 58 targets; residue-contact '
                        'analysis (dock_and_get_interacting_residues) needs a prepared '
                        'receptor PDB on disk: HMGCR, ADRB1, ADRB2, MAOB, DRD2.')
    p.add_argument('--model', default='deepseek-v4-pro',
                   help='Ollama main model name (default: deepseek-v4-pro). '
                        'See --list-models for the built-in options; use the API '
                        'name without the -cloud suffix.')
    think_grp = p.add_mutually_exclusive_group()
    think_grp.add_argument('--think', dest='think', action='store_true', default=None,
                           help="Force Ollama think mode on (overrides auto).")
    think_grp.add_argument('--no-think', dest='think', action='store_false', default=None,
                           help="Force Ollama think mode off (overrides auto).")

    p.add_argument('--ollama-host', default='https://ollama.com', help='Ollama host (default: https://ollama.com).')
    p.add_argument('--ollama-key', default=None, help='Ollama bearer token (or env OLLAMA_API_KEY).')

    p.add_argument('--adversary', choices=['openai', 'anthropic'], default='openai',
                   help='Adversary provider (default: openai).')
    p.add_argument('--adversary-model', default=None,
                   help='Adversary model (default: gpt-5.2 for openai, claude-haiku-4-5-20251001 for anthropic).')
    p.add_argument('--openai-key', default=None, help='OpenAI API key (or env OPENAI_API_KEY).')
    p.add_argument('--anthropic-key', default=None, help='Anthropic API key (or env ANTHROPIC_API_KEY).')
    p.add_argument('--anthropic-base-url', default=None,
                   help='Anthropic API base URL (default: https://api.anthropic.com). '
                        'Given explicitly so an ambient ANTHROPIC_BASE_URL (e.g. a local '
                        'proxy) is bypassed; set this to point at your own proxy if needed.')

    p.add_argument('--context-file', default=os.path.join(_HERE, 'code', 'adversarial_set.md'),
                   help='Initial molecule/score list (default: code/adversarial_set.md). '
                        'Ignored when --resume is given.')
    p.add_argument('--resume', default=None, metavar='JSON',
                   help='Resume a prior run from its JSON messages sidecar (written next '
                        'to the results .md). Skips context-file seeding and the initial '
                        'model turn; continues into the adversary refinement loop from '
                        'the last assistant message already in the sidecar.')
    p.add_argument('--results-dir', default=os.path.join(_HERE, 'results'),
                   help='Where to write the timestamped results .md (default: ./results).')
    p.add_argument('--max-turns', type=int, default=20, help='Safety cap on adversary<->main turns (default: 20).')
    p.add_argument('--max-tool-calls', type=int, default=12,
                   help='Max tool-calling rounds per main-model turn (default: 12). '
                        'Prevents a stuck model from looping forever; once hit, the '
                        'model is forced to emit a text response.')
    p.add_argument('--max-retries', type=int, default=3, help='Retries on Ollama connection errors (default: 3).')
    p.add_argument('--api-timeout', type=float, default=_DEFAULT_API_TIMEOUT,
                   help=f'Per-call API timeout in seconds for the Ollama + adversary clients '
                        f'(default: {_DEFAULT_API_TIMEOUT:.0f}; fail-fast on permanent errors).')
    p.add_argument('--quiet', action='store_true', help='Suppress thinking/content/tool prints.')
    p.add_argument('--rdkit-verbose', action='store_true',
                   help='Re-enable RDKit stderr (SMILES Parse Error) logs. They are silenced '
                        'by default since invalid substituents are now surfaced to the model '
                        'as "invalid SMILES, skipped" entries in the tool results.')
    p.add_argument('--trace', action='store_true',
                   help='Write the full thinking + tool-call trace into the results .md '
                        '(default: only terse section headers + text go to the md; the '
                        'trace otherwise only goes to stdout).')

    p.add_argument('--list-models', action='store_true', help='Print the built-in Ollama model list and exit.')
    p.add_argument('--self-test', action='store_true', help='Run the chemistry tools directly (no LLM keys) and exit.')
    p.add_argument('--mock-tools', action='store_true',
                   help='Replace the three docking-dependent tools (grow_cycle, replace_groups, '
                        'dock_and_get_interacting_residues) with an instant synthetic score, '
                        'skipping real Vina docking. For smoke-testing the LLM/tool-calling and '
                        'adversary wiring in seconds instead of minutes/hours -- not for real runs.')
    return p


def main(argv=None) -> int:
    load_dotenv()
    args = build_arg_parser().parse_args(argv)

    mock_tools.install(args.mock_tools)

    if args.rdkit_verbose:
        RDLogger.EnableLog('rdApp.*')

    if args.list_models:
        print("Ollama models:", ", ".join(OLLAMA_MODELS))
        print("think is auto-disabled for:", ", ".join(sorted(NO_THINK_MODELS)))
        return 0
    if args.self_test:
        self_test()
        return 0

    # Resolve think: explicit flag wins, else auto from model (notebook logic).
    if args.think is None:
        args.think = args.model not in NO_THINK_MODELS

    # Resolve adversary model default per provider.
    if args.adversary_model is None:
        args.adversary_model = 'gpt-5.2' if args.adversary == 'openai' else 'claude-haiku-4-5-20251001'

    if args.resume:
        if not os.path.isfile(args.resume):
            raise SystemExit(f"Resume sidecar not found: {args.resume}")
    elif not os.path.isfile(args.context_file):
        raise SystemExit(f"Context file not found: {args.context_file}")

    run_session(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())