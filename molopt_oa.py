#!/usr/bin/env python3
"""
molopt_oa.py - Command-line adversarial molecule optimization (OpenAI <-> Anthropic).

Variant of molopt.py where the two adversaries are OpenAI and Anthropic instead
of Ollama (main) + OpenAI/Anthropic (critic). `--start` selects which provider
leads: the starter is the tool-calling proposer (it has the chemistry tools),
the other provider is the critique-only adversary. Otherwise the loop is the
same as molopt.py:

  1. The starter reasons over a molecule/docking-score list and may call
     chemistry tools (grow_cycle, replace_groups, make_random_list, related,
     lipinski, dock_and_get_interacting_residues, calculate_SAS_and_NP).
  2. The adversary critiques each proposal.
  3. The starter refines until it replies "Done" (or --max-turns is hit).

Each provider's *native* tool-calling API is used directly (no LangGraph, no
extra deps): OpenAI Chat Completions `tools` and Anthropic Messages `tools`.
The conversation history is kept in the starter's native message format -- the
adversary is a stateless one-shot `critique(prompt) -> str` call, so no
cross-provider message conversion is needed.

Every step is appended to a timestamped Markdown file under --results-dir, and
a JSON messages sidecar (in the starter's native format) is written each turn
for --resume.

Secrets are read from CLI flags or environment variables (optionally a .env
file in the working directory).
"""

import os
import sys
import re
import json
import time
import argparse
import base64

# --- Path + NumPy-2 shim setup (must run before importing the helpers) ------
# The helpers live in ./code; put that dir on the path like molopt.py / test.py do.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, 'code'))

# ODDT 0.7 calls np.in1d, removed in NumPy 2.x. np.isin is a drop-in.
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
# model passes are reported back to it as 'invalid SMILES, skipped' entries in
# the tool results (see MolPropOp.py), so the stderr warnings are redundant.
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Callables the starter model may invoke, in stable order.
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

# Default model names per provider.
OPENAI_DEFAULT_MODEL = 'gpt-5.2'
ANTHROPIC_DEFAULT_MODEL = 'claude-haiku-4-5-20251001'
GEMINI_DEFAULT_MODEL = 'gemini-3-flash-preview'

# --- Prompt text (ported verbatim from molopt.py / notebook cell 11) ---------

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
    """Build the starter model's system message (notebook cell 11 `sys_message`).

    Identical to molopt.py: uses the imported task_specific_prompt /
    task_specific_tools templates.
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


# --- Tool schemas (provider-agnostic JSON Schema; converted per provider) ----

# Each entry: (name, description, parameters JSON Schema). Only the SMILES / count
# argument the tool centers on is marked required; the rest have python defaults
# (predefined substituent sets etc.) so the model can omit them.
_TOOL_DEFS = [
    ('grow_cycle',
     'Add substituents to free carbons of a molecule, dock the results, and return a list of '
     '(molecule SMILES, docking score) pairs. Use it to explore modifications to a promising '
     'molecule. substituents may be a list of SMILES fragments or one of the predefined sets: '
     'e_withdraw, e_donate, withdraw_with_linkers, donate_with_linkers.',
     {'type': 'object',
      'properties': {
          'best_smiles': {'type': 'string', 'description': 'SMILES of the molecule to grow substituents onto.'},
          'best_score': {'type': 'number', 'description': 'Known/estimated docking score of best_smiles (optional).'},
          'substituents': {'type': 'array', 'items': {'type': 'string'},
                           'description': 'Substituent SMILES fragments to add, or a predefined set name (optional).'},
      },
      'required': ['best_smiles']}),
    ('replace_groups',
     'Replace existing substituents in a molecule with new ones, dock the results, and return a '
     'list of (molecule SMILES, docking score) pairs. Use it to test hypotheses about swapping '
     'specific groups. substituents_to_replace / new_substituents may be lists of SMILES fragments '
     'or predefined set names: e_withdraw, e_donate, withdraw_with_linkers, donate_with_linkers.',
     {'type': 'object',
      'properties': {
          'orig_smiles': {'type': 'string', 'description': 'SMILES of the molecule to modify.'},
          'best_score': {'type': 'number', 'description': 'Known/estimated docking score of orig_smiles (optional).'},
          'substituents_to_replace': {'type': 'array', 'items': {'type': 'string'},
                                      'description': 'Groups to replace, or a predefined set name (optional).'},
          'new_substituents': {'type': 'array', 'items': {'type': 'string'},
                               'description': 'Replacement groups, or a predefined set name (optional).'},
      },
      'required': ['orig_smiles']}),
    ('make_random_list',
     'Generate a random list of substituent SMILES fragments of the requested length, drawn from '
     'the predefined electron-withdrawing / electron-donating sets. Useful as input to grow_cycle '
     'or replace_groups when you want a broad sweep of the chemical space.',
     {'type': 'object',
      'properties': {
          'num_items': {'type': 'integer', 'description': 'Number of substituents to select.'},
      },
      'required': ['num_items']}),
    ('related',
     'Generate a list of molecules structurally related to the given molecule(s), returning each '
     'with a few properties. Useful for exploring the chemical space around a promising molecule.',
     {'type': 'object',
      'properties': {
          'smiles_list': {'type': 'array', 'items': {'type': 'string'},
                          'description': 'List of molecule SMILES to find related molecules for.'},
      },
      'required': ['smiles_list']}),
    ('lipinski',
     'Evaluate a list of molecules for drug-likeness via Lipinski\'s rule of five and QED '
     '(0-1, higher = more drug-like). Use it to ensure proposed molecules have good drug-like '
     'properties alongside good docking scores.',
     {'type': 'object',
      'properties': {
          'smiles_list': {'type': 'array', 'items': {'type': 'string'},
                          'description': 'List of molecule SMILES to evaluate.'},
      },
      'required': ['smiles_list']}),
    ('dock_and_get_interacting_residues',
     'Dock a single molecule and return its docking score plus the protein residues it interacts '
     'with and the interaction types. Use only on a molecule already deemed to have a low docking '
     'score and good Lipinski properties, to confirm it binds the expected site.',
     {'type': 'object',
      'properties': {
          'smiles': {'type': 'string', 'description': 'SMILES of the molecule to dock.'},
      },
      'required': ['smiles']}),
    ('calculate_SAS_and_NP',
     'Calculate the synthetic accessibility score (SAS, 1 easy - 10 hard) and natural-product '
     'likeness score (NP, -5 to 5, higher = more natural-product-like) for a list of molecules. '
     'Call for promising molecules with good docking scores and good Lipinski properties.',
     {'type': 'object',
      'properties': {
          'smiles_list': {'type': 'array', 'items': {'type': 'string'},
                          'description': 'List of molecule SMILES to score.'},
      },
      'required': ['smiles_list']}),
]


def _openai_tools():
    """OpenAI Chat Completions tool spec list."""
    return [{'type': 'function', 'function': {'name': n, 'description': d, 'parameters': p}}
            for (n, d, p) in _TOOL_DEFS]


def _anthropic_tools():
    """Anthropic Messages tool spec list."""
    return [{'name': n, 'description': d, 'input_schema': p} for (n, d, p) in _TOOL_DEFS]


# --- Gemini tool schema (google-genai FunctionDeclarations) ------------------
#
# Gemini's schema uses uppercase type names (OBJECT/STRING/NUMBER/...) instead
# of JSON Schema's lowercase ones, so we recurse over the provider-agnostic
# _TOOL_DEFS and lift each sub-schema into a types.Schema.

_GEMINI_TYPE_MAP = {
    'object': 'OBJECT', 'string': 'STRING', 'number': 'NUMBER',
    'integer': 'INTEGER', 'boolean': 'BOOLEAN', 'array': 'ARRAY',
}


def _schema_to_gemini(schema):
    """Recursively convert a JSON Schema dict into a google.genai types.Schema."""
    from google.genai import types
    if not isinstance(schema, dict):
        return schema
    kwargs = {}
    t = schema.get('type')
    if t:
        kwargs['type'] = _GEMINI_TYPE_MAP.get(t, str(t).upper())
    if 'description' in schema:
        kwargs['description'] = schema['description']
    if 'properties' in schema:
        kwargs['properties'] = {k: _schema_to_gemini(v)
                                for k, v in schema['properties'].items()}
    if 'required' in schema:
        kwargs['required'] = list(schema['required'])
    if 'items' in schema:
        kwargs['items'] = _schema_to_gemini(schema['items'])
    return types.Schema(**kwargs)


def _gemini_tools():
    """A single google.genai Tool holding all function declarations."""
    from google.genai import types
    decls = [types.FunctionDeclaration(name=n, description=d,
                                       parameters=_schema_to_gemini(p))
             for (n, d, p) in _TOOL_DEFS]
    return [types.Tool(function_declarations=decls)]


# --- Minimal .env loader (no extra dependency) ------------------------------

_ENV_KEYS = ('OPENAI_API_KEY', 'ANTHROPIC_API_KEY', 'GEMINI_API_KEY')


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


# --- Shared helpers ----------------------------------------------------------

def _vprint(verbose: bool, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)


_DONE_RE = re.compile(r'Done[.!?]*\s*$', re.IGNORECASE)


def is_done(text: str) -> bool:
    """True if the model signalled completion with "Done" (ported from molopt.py)."""
    if not text or not text.strip():
        return False
    s = text.strip()
    if _DONE_RE.fullmatch(s):
        return True
    last_line = s.splitlines()[-1].strip()
    return bool(_DONE_RE.fullmatch(last_line))


def _run_tool(fn_name: str, fn_args: dict, verbose: bool):
    """Execute a tool call, returning a string result (errors are fed back, not raised)."""
    if fn_name not in AVAILABLE_FUNCTIONS:
        _vprint(verbose, f"  [unknown tool {fn_name}; ignored]")
        return f"Error: unknown tool '{fn_name}'."
    try:
        result = AVAILABLE_FUNCTIONS[fn_name](**fn_args)
    except Exception as err:
        result = (f"Tool '{fn_name}' raised an error: {err}. "
                  f"Check the argument names/types and retry.")
        _vprint(verbose, f"  [tool error: {err}]")
    _vprint(verbose, f"Result: {result}")
    _vprint(verbose, '-' * 72)
    return result


# --- API-error classification ------------------------------------------------
#
# Both provider SDKs default to max_retries=2 with exponential backoff and a
# 600s per-call timeout. That is fine for a transient 5xx, but for a PERMANENT
# error -- bad auth (401/403) or an out-of-credits 429 -- the SDK's backoff
# (honouring Retry-After) makes the call appear to hang for many minutes with
# 0% CPU and no output, instead of failing fast. We therefore build both
# clients with max_retries=0 + an explicit timeout (the manual retry loops in
# chat_turn control retries), and refuse to retry errors that can never succeed.

_DEFAULT_API_TIMEOUT = 120.0  # seconds per call; well under the SDK's 600s default


def _is_retryable(err) -> bool:
    """True for transient errors worth retrying; False for permanent ones.

    Permanent (do not retry): any 4xx client error (400 bad request / invalid
    API key, 401/403 auth), and a 429 that is actually billing/quota exhaustion
    rather than a transient rate limit. Connection errors and 5xx have no
    status code and are treated as retryable.

    Provider-agnostic: OpenAI/Anthropic surface ``status_code``, while
    google-genai raises ``ClientError`` carrying ``code`` (a bad Gemini key is a
    400 INVALID_ARGUMENT, not a 401). We check both attributes.
    """
    sc = getattr(err, 'status_code', None)
    if sc is None:
        sc = getattr(err, 'code', None)
    if sc is not None and 400 <= sc < 500:
        if sc == 429:
            # Transient rate limit vs. permanent quota/billing exhaustion.
            low = str(err).lower()
            if any(k in low for k in ('credit', 'quota', 'billing', 'insufficient')):
                return False
            return True
        return False  # 4xx (bad request / auth): retrying cannot help
    return True  # 5xx / connection error / unknown: retryable


# --- OpenAI actor (Chat Completions, native tool calling) -------------------

class OpenAIActor:
    """OpenAI actor. Messages are kept in Chat Completions format (plain dicts):
    {'role':'system'|'user'|'assistant'|'tool', 'content':..., 'tool_calls':...,
     'tool_call_id':...}. The system message is the first entry.

    Implements both `chat_turn` (tool-calling proposer) and `critique`
    (stateless adversary), so the same class can serve as starter or adversary.
    """

    def __init__(self, model: str, api_key: str, timeout: float = _DEFAULT_API_TIMEOUT):
        from openai import OpenAI  # lazy import
        self.model = model
        # max_retries=0: the manual loop in chat_turn handles retries, and we do
        # NOT want the SDK's internal backoff (it hangs for minutes on a 429).
        self.client = OpenAI(api_key=api_key, timeout=timeout, max_retries=0)
        self.tools = _openai_tools()

    # -- adversary (no tools, single shot) --
    def critique(self, prompt: str) -> str:
        resp = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {'role': 'system', 'content': ADVERSARY_INSTRUCTIONS},
                {'role': 'user', 'content': prompt},
            ],
        )
        return (resp.choices[0].message.content or '').strip()

    # -- proposer (tool-calling loop) --
    def chat_turn(self, messages, prompt, *, max_retries, max_tool_calls, verbose):
        """Append a user prompt and run the model until it stops calling tools.

        Returns (messages, last_assistant_text, trace). Mirrors molopt.py's
        chat_turn: connection-error retry, tool-error isolation, per-turn
        tool-call cap that forces a text-only summary once hit.
        """
        messages = list(messages)
        messages.append({'role': 'user', 'content': prompt})

        trace = []
        tool_rounds = 0
        while True:
            force_tools = self.tools if tool_rounds < max_tool_calls else []
            if tool_rounds >= max_tool_calls:
                _vprint(verbose, f"  [tool-call cap ({max_tool_calls}) reached; "
                                 f"forcing a text-only response]")

            response = None
            last_err = None
            for attempt in range(1, max_retries + 1):
                try:
                    response = self.client.chat.completions.create(
                        model=self.model,
                        messages=messages,
                        tools=force_tools,
                    )
                    break
                except Exception as err:
                    last_err = err
                    _vprint(verbose, f"  [openai error, attempt {attempt}/{max_retries}: {err}]")
                    if not _is_retryable(err):
                        raise  # permanent (auth / out-of-credits): fail fast, don't loop
                    if attempt < max_retries:
                        time.sleep(2 * attempt)
            if response is None:
                msg = (f"The previous call failed with a connection error: {last_err}. "
                       f"Please proceed from the last step.")
                _vprint(verbose, f"  [giving up after {max_retries} retries; asking model to continue]")
                messages.append({'role': 'user', 'content': msg})
                continue

            choice = response.choices[0]
            assistant_msg = choice.message
            # Build a JSON-safe dict for the sidecar. tool_calls carry their
            # server-assigned ids, which the matching tool messages must reference.
            entry = {'role': 'assistant', 'content': assistant_msg.content or ''}
            tool_calls = []
            for tc in (assistant_msg.tool_calls or []):
                try:
                    args = json.loads(tc.function.arguments) if tc.function.arguments else {}
                except Exception:
                    args = {}
                tool_calls.append({'id': tc.id, 'name': tc.function.name, 'arguments': args})
                _vprint(verbose, f"Calling {tc.function.name} with arguments {args}")
            if tool_calls:
                entry['tool_calls'] = [{'id': tc['id'], 'type': 'function',
                                        'function': {'name': tc['name'],
                                                     'arguments': json.dumps(tc['arguments'])}}
                                       for tc in tool_calls]
            messages.append(entry)
            _vprint(verbose, '-' * 72)
            _vprint(verbose, "Content: ", entry['content'])
            _vprint(verbose, '-' * 72)
            trace.append(f"### content\n{entry['content']}")

            if tool_calls and tool_rounds < max_tool_calls:
                tool_rounds += 1
                for tc in tool_calls:
                    result = _run_tool(tc['name'], tc['arguments'], verbose)
                    trace.append(f"### tool call: {tc['name']}\nargs: {tc['arguments']}\nresult: {result}")
                    messages.append({'role': 'tool', 'tool_call_id': tc['id'], 'content': str(result)})
            else:
                if tool_calls:
                    # Over cap. OpenAI requires every assistant tool_call to be
                    # answered by a matching `tool` message before the next call,
                    # so emit stub results for the calls we are NOT executing,
                    # then nudge for a text-only summary with no tools available.
                    _vprint(verbose, "  [over cap; requesting a text summary without tools]")
                    for tc in tool_calls:
                        messages.append({'role': 'tool', 'tool_call_id': tc['id'],
                                         'content': 'Tool-call limit reached for this turn; '
                                                    'this call was not executed.'})
                    messages.append({'role': 'user', 'content':
                        "You have reached the tool-call limit for this turn. Do not call any "
                        "more tools. Summarize your best proposed molecules so far, their "
                        "estimated docking scores, and your reasoning. Do NOT say 'Done' — "
                        "you will receive adversary feedback next and then refine your proposals."})
                    try:
                        final = self.client.chat.completions.create(
                            model=self.model, messages=messages, tools=[])
                        fm = final.choices[0].message
                        messages.append({'role': 'assistant', 'content': fm.content or ''})
                        _vprint(verbose, '-' * 72)
                        _vprint(verbose, "Content: ", fm.content)
                        _vprint(verbose, '-' * 72)
                        trace.append(f"### content (over-cap summary)\n{fm.content or ''}")
                    except Exception as err:
                        _vprint(verbose, f"  [final summary call failed: {err}]")
                break  # turn is done

        return messages, self.last_assistant_text(messages), trace

    @staticmethod
    def last_assistant_text(messages) -> str:
        """Most recent non-empty assistant message text (content can be empty
        when the assistant message only carried tool_calls)."""
        for m in reversed(messages):
            if m.get('role') == 'assistant':
                content = m.get('content') or ''
                if content and content.strip():
                    return content
        return ''


# --- Anthropic actor (Messages API, native tool calling) --------------------

class AnthropicActor:
    """Anthropic actor. Messages are kept in Messages-API format (plain dicts):
    {'role':'system'|'user'|'assistant', 'content': str | list[blocks]}. The
    system message is the first entry and is passed via the `system=` parameter
    on each call (not inside `messages`). Tool results are sent as user-role
    `tool_result` content blocks.

    Implements both `chat_turn` (tool-calling proposer) and `critique`
    (stateless adversary), so the same class can serve as starter or adversary.
    """

    def __init__(self, model: str, api_key: str, base_url: str = "https://api.anthropic.com",
                 timeout: float = _DEFAULT_API_TIMEOUT):
        from anthropic import Anthropic  # lazy import
        self.model = model
        # Explicit base_url: bypass any ambient ANTHROPIC_BASE_URL (e.g. a local
        # proxy that doesn't serve Anthropic model names), routing to the real API.
        # max_retries=0 + timeout: see _is_retryable -- avoid the SDK's minutes-long
        # backoff on a permanent 429/401; chat_turn's manual loop handles retries.
        self.client = Anthropic(api_key=api_key, base_url=base_url,
                                timeout=timeout, max_retries=0)
        self.tools = _anthropic_tools()

    def _split_system(self, messages):
        """Pull a leading system message out of the list -> (system_str, rest)."""
        if messages and messages[0].get('role') == 'system':
            sysmsg = messages[0].get('content') or ''
            if not isinstance(sysmsg, str):
                # Flatten a block list to text just in case.
                sysmsg = ''.join(b.get('text', '') for b in sysmsg if isinstance(b, dict))
            return sysmsg, list(messages[1:])
        return build_system_message(), list(messages)

    # -- adversary (no tools, single shot) --
    def critique(self, prompt: str) -> str:
        resp = self.client.messages.create(
            model=self.model,
            system=ADVERSARY_INSTRUCTIONS,
            max_tokens=4096,
            messages=[{"role": "user", "content": prompt}],
        )
        parts = [b.text for b in resp.content if getattr(b, "type", None) == "text"]
        return "".join(parts).strip()

    # -- proposer (tool-calling loop) --
    def chat_turn(self, messages, prompt, *, max_retries, max_tool_calls, verbose):
        """Append a user prompt and run the model until it stops calling tools.

        Returns (messages, last_assistant_text, trace). Mirrors molopt.py's
        chat_turn: connection-error retry, tool-error isolation, per-turn
        tool-call cap that forces a text-only summary once hit.
        """
        messages = list(messages)
        messages.append({'role': 'user', 'content': prompt})

        trace = []
        tool_rounds = 0
        while True:
            system, api_messages = self._split_system(messages)
            force_tools = self.tools if tool_rounds < max_tool_calls else []
            if tool_rounds >= max_tool_calls:
                _vprint(verbose, f"  [tool-call cap ({max_tool_calls}) reached; "
                                 f"forcing a text-only response]")

            response = None
            last_err = None
            for attempt in range(1, max_retries + 1):
                try:
                    response = self.client.messages.create(
                        model=self.model,
                        system=system,
                        max_tokens=4096,
                        messages=api_messages,
                        tools=force_tools,
                    )
                    break
                except Exception as err:
                    last_err = err
                    _vprint(verbose, f"  [anthropic error, attempt {attempt}/{max_retries}: {err}]")
                    if not _is_retryable(err):
                        raise  # permanent (auth / out-of-credits): fail fast, don't loop
                    if attempt < max_retries:
                        time.sleep(2 * attempt)
            if response is None:
                msg = (f"The previous call failed with a connection error: {last_err}. "
                       f"Please proceed from the last step.")
                _vprint(verbose, f"  [giving up after {max_retries} retries; asking model to continue]")
                messages.append({'role': 'user', 'content': msg})
                continue

            # Rebuild a JSON-safe assistant entry from the content blocks.
            blocks = []
            text_parts = []
            tool_uses = []
            for b in response.content:
                btype = getattr(b, 'type', None)
                if btype == 'text':
                    blocks.append({'type': 'text', 'text': b.text})
                    text_parts.append(b.text)
                elif btype == 'tool_use':
                    blocks.append({'type': 'tool_use', 'id': b.id, 'name': b.name, 'input': b.input})
                    tool_uses.append({'id': b.id, 'name': b.name, 'input': b.input})
                    _vprint(verbose, f"Calling {b.name} with arguments {b.input}")
            entry = {'role': 'assistant', 'content': blocks if blocks else ''}
            messages.append(entry)
            content_str = ''.join(text_parts)
            _vprint(verbose, '-' * 72)
            _vprint(verbose, "Content: ", content_str)
            _vprint(verbose, '-' * 72)
            trace.append(f"### content\n{content_str}")

            if tool_uses and tool_rounds < max_tool_calls:
                tool_rounds += 1
                tool_results = []
                for tu in tool_uses:
                    result = _run_tool(tu['name'], tu['input'], verbose)
                    trace.append(f"### tool call: {tu['name']}\nargs: {tu['input']}\nresult: {result}")
                    tool_results.append({'type': 'tool_result', 'tool_use_id': tu['id'], 'content': str(result)})
                messages.append({'role': 'user', 'content': tool_results})
            else:
                if tool_uses:
                    # Over cap. Anthropic requires every assistant tool_use block
                    # to be answered by a matching tool_result in the next user
                    # message, and messages must strictly alternate user/assistant
                    # -- so we emit the stub results AND the nudge text together in
                    # one user message, then call with no tools for a text summary.
                    _vprint(verbose, "  [over cap; requesting a text summary without tools]")
                    nudge = ("You have reached the tool-call limit for this turn. Do not call any "
                             "more tools. Summarize your best proposed molecules so far, their "
                             "estimated docking scores, and your reasoning. Do NOT say 'Done' — "
                             "you will receive adversary feedback next and then refine your proposals.")
                    over_cap_blocks = [{'type': 'tool_result', 'tool_use_id': tu['id'],
                                        'content': 'Tool-call limit reached for this turn; '
                                                   'this call was not executed.',
                                        'is_error': False} for tu in tool_uses]
                    over_cap_blocks.append({'type': 'text', 'text': nudge})
                    messages.append({'role': 'user', 'content': over_cap_blocks})
                    try:
                        system, api_messages = self._split_system(messages)
                        final = self.client.messages.create(
                            model=self.model, system=system, max_tokens=4096,
                            messages=api_messages, tools=[])
                        fblocks = []
                        ftext = []
                        for b in final.content:
                            if getattr(b, 'type', None) == 'text':
                                fblocks.append({'type': 'text', 'text': b.text})
                                ftext.append(b.text)
                        messages.append({'role': 'assistant', 'content': fblocks if fblocks else ''})
                        _vprint(verbose, '-' * 72)
                        _vprint(verbose, "Content: ", ''.join(ftext))
                        _vprint(verbose, '-' * 72)
                        trace.append(f"### content (over-cap summary)\n{''.join(ftext)}")
                    except Exception as err:
                        _vprint(verbose, f"  [final summary call failed: {err}]")
                break  # turn is done

        return messages, self.last_assistant_text(messages), trace

    @staticmethod
    def last_assistant_text(messages) -> str:
        """Most recent non-empty assistant text. Anthropic assistant content is
        a list of blocks; pull the text blocks. (A tool_use-only assistant
        message has no text -> scan back for the latest one that does.)"""
        for m in reversed(messages):
            if m.get('role') != 'assistant':
                continue
            content = m.get('content')
            if isinstance(content, str):
                if content.strip():
                    return content
            elif isinstance(content, list):
                text = ''.join(b.get('text', '') for b in content
                               if isinstance(b, dict) and b.get('type') == 'text')
                if text.strip():
                    return text
        return ''


# --- Gemini actor (google-genai, native function calling) --------------------

# The google-genai SDK prints a benign "Direct use of automatic function
# calling ... is not recommended" notice on every call. We do MANUAL function
# calling (FunctionDeclarations + our own function_response handling), which
# the notice does not apply to -- silence it so it does not spam the log.
import warnings as _warnings
_warnings.filterwarnings('ignore', message='.*automatic function calling.*')


def _ts_encode(blob) -> str | None:
    """thought_signature is opaque bytes (len ~hundreds). Base64-encode it so it
    survives json.dump into the sidecar (the sidecar uses default=str, which
    would stringify bytes irreversibly). Returns None when there is no
    signature (non-thinking Gemini models, or a part without one)."""
    if blob is None:
        return None
    if isinstance(blob, str):
        blob = blob.encode('utf-8', 'surrogateescape')
    return base64.b64encode(blob).decode('ascii')


def _ts_decode(sig) -> bytes | None:
    """Inverse of _ts_encode -> bytes for the SDK Part.thought_signature field."""
    if not sig:
        return None
    return base64.b64decode(sig)


class GeminiActor:
    """Google Gemini actor (google-genai). Messages are kept in a JSON-safe
    Gemini-shaped dict format::

        {'role': 'system'|'user'|'model', 'content': str}   # system only
        {'role': 'user'|'model', 'parts': [
            {'text': ...} |
            {'function_call': {'name': ..., 'args': {...}}} |
            {'function_response': {'name': ..., 'response': {...}}}]}

    The system message is the first entry (passed via system_instruction on
    each call, like Anthropic). Tool results are sent back as user-role
    ``function_response`` parts. Implements both ``chat_turn`` (tool-calling
    proposer) and ``critique`` (stateless adversary).
    """

    def __init__(self, model: str, api_key: str, timeout: float = _DEFAULT_API_TIMEOUT):
        from google import genai
        from google.genai import types
        self.model = model
        # HttpOptions.timeout is in milliseconds (the SDK rejects a deadline
        # below 10s). The manual loop + _is_retryable handle retries / fail-fast
        # on permanent 4xx (bad key, bad request).
        self.client = genai.Client(
            api_key=api_key,
            http_options=types.HttpOptions(timeout=float(timeout) * 1000))
        self.tools = _gemini_tools()

    # -- message-format helpers --
    def _split_system(self, messages):
        """Pull a leading system message out -> (system_str, rest)."""
        if messages and messages[0].get('role') == 'system':
            sysmsg = messages[0].get('content') or ''
            if not isinstance(sysmsg, str):
                sysmsg = ''.join(b.get('text', '') for b in sysmsg
                                 if isinstance(b, dict))
            return sysmsg, list(messages[1:])
        return build_system_message(), list(messages)

    def _to_contents(self, msg_dicts):
        """Convert the JSON-safe message dicts into SDK Content objects.

        Gemini 3.x thinking models attach a ``thought_signature`` (opaque bytes)
        to function_call parts (and sometimes the final response part). It MUST
        be echoed back with the matching function_response, and the model's own
        function_call turn must carry it when resent -- otherwise the API
        rejects the request with 400 "Function call is missing a
        thought_signature". We store it base64-encoded at the part-dict level
        and decode it back to bytes here.
        """
        from google.genai import types
        out = []
        for m in msg_dicts:
            role = m['role']
            parts = []
            for p in m.get('parts', []):
                ts = _ts_decode(p.get('thought_signature'))
                if 'text' in p:
                    parts.append(types.Part(text=p['text'], thought_signature=ts))
                elif 'function_call' in p:
                    fc = p['function_call']
                    parts.append(types.Part(function_call=types.FunctionCall(
                        name=fc['name'], args=fc.get('args') or {}),
                        thought_signature=ts))
                elif 'function_response' in p:
                    fr = p['function_response']
                    parts.append(types.Part(function_response=types.FunctionResponse(
                        name=fr['name'], response=fr.get('response') or {}),
                        thought_signature=ts))
            out.append(types.Content(role=role, parts=parts))
        return out

    # -- adversary (no tools, single shot) --
    def critique(self, prompt: str) -> str:
        from google.genai import types
        cfg = types.GenerateContentConfig(system_instruction=ADVERSARY_INSTRUCTIONS)
        r = self.client.models.generate_content(
            model=self.model,
            contents=[types.Content(role='user', parts=[types.Part(text=prompt)])],
            config=cfg)
        return ''.join(getattr(p, 'text', '') or ''
                       for p in r.candidates[0].content.parts).strip()

    # -- proposer (tool-calling loop) --
    def chat_turn(self, messages, prompt, *, max_retries, max_tool_calls, verbose):
        """Append a user prompt and run the model until it stops calling tools.

        Returns (messages, last_assistant_text, trace). Mirrors the other actors:
        connection-error retry, tool-error isolation, per-turn tool-call cap
        that forces a text-only summary once hit.
        """
        from google.genai import types
        messages = list(messages)
        messages.append({'role': 'user', 'parts': [{'text': prompt}]})

        trace = []
        tool_rounds = 0
        while True:
            system, rest = self._split_system(messages)
            force_tools = self.tools if tool_rounds < max_tool_calls else []
            if tool_rounds >= max_tool_calls:
                _vprint(verbose, f"  [tool-call cap ({max_tool_calls}) reached; "
                                 f"forcing a text-only response]")

            cfg_kwargs = {'system_instruction': system}
            if force_tools:
                cfg_kwargs['tools'] = force_tools
            cfg = types.GenerateContentConfig(**cfg_kwargs)

            response = None
            last_err = None
            for attempt in range(1, max_retries + 1):
                try:
                    response = self.client.models.generate_content(
                        model=self.model, contents=self._to_contents(rest), config=cfg)
                    break
                except Exception as err:
                    last_err = err
                    _vprint(verbose, f"  [gemini error, attempt {attempt}/{max_retries}: {err}]")
                    if not _is_retryable(err):
                        raise  # permanent (bad key / bad request): fail fast
                    if attempt < max_retries:
                        time.sleep(2 * attempt)
            if response is None:
                msg = (f"The previous call failed with a connection error: {last_err}. "
                       f"Please proceed from the last step.")
                _vprint(verbose, f"  [giving up after {max_retries} retries; asking model to continue]")
                messages.append({'role': 'user', 'parts': [{'text': msg}]})
                continue

            parts = response.candidates[0].content.parts
            sdk_parts = []
            text_parts = []
            tool_calls = []
            for p in parts:
                # Capture the Part's thought_signature (Gemini 3.x thinking models).
                # It is opaque bytes that MUST be echoed back on the matching
                # function_response and preserved on the resent model turn, so
                # we base64-encode it into the part dict (see _to_contents).
                sig = _ts_encode(getattr(p, 'thought_signature', None))
                fc = getattr(p, 'function_call', None)
                if fc:
                    args = dict(fc.args) if fc.args else {}
                    tool_calls.append({'name': fc.name, 'args': args, 'sig': sig})
                    sdk_parts.append({'function_call': {'name': fc.name, 'args': args},
                                      'thought_signature': sig})
                    _vprint(verbose, f"Calling {fc.name} with arguments {args}")
                else:
                    txt = getattr(p, 'text', None) or ''
                    if txt:
                        text_parts.append(txt)
                        sdk_parts.append({'text': txt, 'thought_signature': sig})
            messages.append({'role': 'model', 'parts': sdk_parts})
            content_str = ''.join(text_parts)
            _vprint(verbose, '-' * 72)
            _vprint(verbose, "Content: ", content_str)
            _vprint(verbose, '-' * 72)
            trace.append(f"### content\n{content_str}")

            if tool_calls and tool_rounds < max_tool_calls:
                tool_rounds += 1
                fr_parts = []
                for tc in tool_calls:
                    result = _run_tool(tc['name'], tc['args'], verbose)
                    trace.append(f"### tool call: {tc['name']}\nargs: {tc['args']}\nresult: {result}")
                    # Echo the function_call's thought_signature on the response.
                    fr_parts.append({'function_response': {'name': tc['name'],
                                                            'response': {'result': str(result)}},
                                      'thought_signature': tc.get('sig')})
                messages.append({'role': 'user', 'parts': fr_parts})
            else:
                if tool_calls:
                    # Over cap. Gemini requires every function_call to be answered
                    # by a matching function_response in the next user turn, so emit
                    # stub responses for the calls we are NOT executing plus a nudge
                    # for a text-only summary, then call with no tools available.
                    _vprint(verbose, "  [over cap; requesting a text summary without tools]")
                    nudge = ("You have reached the tool-call limit for this turn. Do not call any "
                             "more tools. Summarize your best proposed molecules so far, their "
                             "estimated docking scores, and your reasoning. Do NOT say 'Done' — "
                             "you will receive adversary feedback next and then refine your proposals.")
                    over_cap = [{'function_response': {'name': tc['name'],
                                                       'response': {'result':
                           'Tool-call limit reached for this turn; this call was not executed.'}},
                                 'thought_signature': tc.get('sig')}
                                for tc in tool_calls]
                    over_cap.append({'text': nudge})
                    messages.append({'role': 'user', 'parts': over_cap})
                    try:
                        system, rest = self._split_system(messages)
                        final = self.client.models.generate_content(
                            model=self.model, contents=self._to_contents(rest),
                            config=types.GenerateContentConfig(system_instruction=system))
                        fparts = []
                        ftext = []
                        for p in final.candidates[0].content.parts:
                            txt = getattr(p, 'text', None) or ''
                            if txt:
                                fparts.append({'text': txt})
                                ftext.append(txt)
                        messages.append({'role': 'model', 'parts': fparts})
                        _vprint(verbose, '-' * 72)
                        _vprint(verbose, "Content: ", ''.join(ftext))
                        _vprint(verbose, '-' * 72)
                        trace.append(f"### content (over-cap summary)\n{''.join(ftext)}")
                    except Exception as err:
                        _vprint(verbose, f"  [final summary call failed: {err}]")
                break  # turn is done

        return messages, self.last_assistant_text(messages), trace

    @staticmethod
    def last_assistant_text(messages) -> str:
        """Most recent non-empty model text. A function_call-only model message
        has no text -> scan back for the latest one that does."""
        for m in reversed(messages):
            if m.get('role') != 'model':
                continue
            for p in m.get('parts', []):
                if 'text' in p and p['text'].strip():
                    return p['text']
        return ''


# --- Actor construction ------------------------------------------------------

def _build_actor(provider, model, key, args):
    """Construct the actor for a provider name."""
    if provider == 'openai':
        return OpenAIActor(model, key, timeout=args.api_timeout)
    if provider == 'anthropic':
        base_url = args.anthropic_base_url or "https://api.anthropic.com"
        return AnthropicActor(model, key, base_url=base_url, timeout=args.api_timeout)
    if provider == 'gemini':
        return GeminiActor(model, key, timeout=args.api_timeout)
    raise ValueError(f"unknown provider: {provider}")


def make_actors(args):
    """Build (starter, adversary, starter_label, adversary_label) from args.

    starter: the tool-calling proposer (OpenAIActor / AnthropicActor / GeminiActor).
    adversary: the critique-only actor of a *different* provider, chosen by
    --adversary (defaulting to the natural counterpart: openai<->anthropic pair
    up, gemini defaults to an openai adversary to match the Ollama study).
    Only the two providers actually used need their API keys set.
    """
    keys = {
        'openai': args.openai_key or os.environ.get('OPENAI_API_KEY') or '',
        'anthropic': args.anthropic_key or os.environ.get('ANTHROPIC_API_KEY') or '',
        'gemini': args.gemini_key or os.environ.get('GEMINI_API_KEY') or '',
    }
    models = {'openai': args.openai_model, 'anthropic': args.anthropic_model,
              'gemini': args.gemini_model}

    start = args.start
    if args.adversary:
        adversary = args.adversary
    elif start == 'openai':
        adversary = 'anthropic'
    elif start == 'anthropic':
        adversary = 'openai'
    else:  # gemini
        adversary = 'openai'
    if adversary == start and not args.adversary:
        # Only block the *default-derived* same-provider case (would indicate a bug in the
        # fallback table above). An explicit --adversary matching --start is allowed -- it's
        # a same-provider self-critique run (e.g. openai proposer vs. openai adversary), used
        # to remove the cross-provider-adversary confound when comparing against other sets.
        raise SystemExit(f"--adversary ({adversary}) must differ from --start ({start}).")

    if not keys[start]:
        raise SystemExit(f"{start} needs --{start}-key or {start.upper()}_API_KEY.")
    if not keys[adversary]:
        raise SystemExit(f"{adversary} (adversary) needs --{adversary}-key or "
                         f"{adversary.upper()}_API_KEY.")

    starter = _build_actor(start, models[start], keys[start], args)
    adv = _build_actor(adversary, models[adversary], keys[adversary], args)
    label = lambda p: f'{p}/{models[p]}'
    return starter, adv, label(start), label(adversary)


# --- Messages sidecar (for --resume) ----------------------------------------

def _serialize_messages(messages) -> list:
    """The starter's messages are already plain JSON-safe dicts (both providers
    build dict entries, never provider SDK objects), so they serialize as-is."""
    return [dict(m) for m in messages]


def write_sidecar(sidecar_path: str, *, args, messages, written_at_turn, status,
                  starter_label, adversary_label) -> None:
    """Atomically write the messages list + run metadata as a JSON sidecar.

    Written incrementally (after each turn) so a run killed mid-way still
    leaves a usable sidecar of the last completed turn. Atomic via tmp+rename.
    """
    payload = {
        'protein': args.protein,
        'start': args.start,
        'openai_model': args.openai_model,
        'anthropic_model': args.anthropic_model,
        'gemini_model': args.gemini_model,
        'starter': starter_label,
        'adversary': adversary_label,
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


def load_sidecar(path: str, starter) -> list:
    """Load the messages list from a JSON sidecar, ready to feed back to the
    starter. Ensures the first message is a system prompt; if the sidecar is
    missing one (e.g. an older/edited file), prepend the current system message.
    """
    with open(path, 'r') as f:
        payload = json.load(f)
    messages = [dict(m) for m in payload.get('messages', [])]
    if not messages or messages[0].get('role') != 'system':
        messages.insert(0, {'role': 'system', 'content': build_system_message()})
    return messages


# --- Session runner ---------------------------------------------------------

def run_session(args) -> str:
    """Run the full adversarial session. Returns the results file path."""
    # Configure the shared scoring state (same mutable object the helpers use).
    scoring_args[0] = os.cpu_count()
    scoring_args[1] = args.protein

    starter, adversary, starter_label, adversary_label = make_actors(args)

    os.makedirs(args.results_dir, exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    safe = args.start.replace('/', '-')
    results_path = os.path.join(args.results_dir, f"{safe}_{args.protein}_{timestamp}.md")
    # JSON messages sidecar lives next to the results .md (same stem). Always-on;
    # written incrementally so a stopped/killed run can be resumed via --resume.
    sidecar_path = os.path.splitext(results_path)[0] + '.json'
    with open(results_path, 'w') as f:
        f.write(f'# Adversarial Design Session (OpenAI<->Anthropic) - {timestamp}\n')
        f.write(f'# protein: {args.protein} | starter: {starter_label} (tools) '
                f'| adversary: {adversary_label} (critique)\n')
        if args.resume:
            f.write(f'# resumed from: {args.resume}\n')
        f.write('\n')

    def log(section, text):
        with open(results_path, 'a') as f:
            f.write(f'\n{section}\n{text}\n')

    if args.resume:
        # Resume: seed the conversation from a prior run's JSON sidecar instead
        # of the context file, and skip the initial starter turn. We pick up from
        # the last assistant message already in the sidecar.
        print(f"Resuming from sidecar -> {args.resume}")
        messages = load_sidecar(args.resume, starter)
        last = starter.last_assistant_text(messages)
        log('# Resumed from sidecar:', args.resume)
        log('# Last assistant text at resume:', last or '(none)')
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=0, status='resumed',
                      starter_label=starter_label, adversary_label=adversary_label)
        turn = 0
    else:
        messages = [{'role': 'system', 'content': build_system_message()}]

        # Initial prompt: the molecule / docking-score list for this protein.
        with open(args.context_file, 'r') as f:
            context = f.read()
        first_prompt = f'\n  Here is a list of molecules and their docking scores:\n  {context}\n'

        print(f"Starting session -> {results_path}")
        messages, last, trace = starter.chat_turn(
            messages, first_prompt,
            max_retries=args.max_retries, max_tool_calls=args.max_tool_calls,
            verbose=not args.quiet,
        )
        log('# Initial model response:', last)
        if args.trace:
            log('# Trace:', '\n'.join(trace))
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=0, status='in_progress',
                      starter_label=starter_label, adversary_label=adversary_label)
        turn = 0

    while not is_done(last) and turn < args.max_turns:
        turn += 1
        print(f"\n=== Turn {turn}/{args.max_turns} ===")
        if last and last.strip():
            try:
                adv = adversary.critique(last)
                log('# Adversary feedback:', adv)
                _vprint(not args.quiet, f"[adversary {adversary_label} "
                                        f"replied ({len(adv)} chars)]")
            except Exception as err:
                # A permanent error (bad auth / out-of-credits) will recur every
                # turn, so a "rescued" run would be silently critique-free -- fail
                # fast instead. A transient error (rate limit, 5xx, connection) is
                # recovered in-line so one blip doesn't kill a long run.
                if not _is_retryable(err):
                    raise
                adv = (f"The adversary model could not be reached (error: {err}). "
                       f"Review your latest proposal yourself, correct any flaws you "
                       f"can identify, and present your best molecules with estimated "
                       f"scores. Say 'Done' if you are finished.")
                log('# Adversary feedback: [unavailable]', adv)
                _vprint(not args.quiet, f"[adversary error: {err}]")
        else:
            # Starter ended its turn with no text; nudge it to summarize instead
            # of calling the adversary with empty input.
            adv = ("Your last response had no text. Summarize your best proposed "
                   "molecules, their estimated docking scores, and your reasoning. "
                   "Say 'Done' if you are finished.")
            log('# Adversary feedback: [skipped - empty starter response]', adv)
            _vprint(not args.quiet, "[adversary skipped - empty starter response]")

        messages, last, turn_trace = starter.chat_turn(
            messages, adv,
            max_retries=args.max_retries, max_tool_calls=args.max_tool_calls,
            verbose=not args.quiet,
        )
        log('# Model response:', last)
        if args.trace:
            log('# Trace:', '\n'.join(turn_trace))
        write_sidecar(sidecar_path, args=args, messages=messages,
                      written_at_turn=turn, status='in_progress',
                      starter_label=starter_label, adversary_label=adversary_label)

    status = 'Done' if is_done(last) else f'MAX_TURNS_REACHED (last={last!r})'
    with open(results_path, 'a') as f:
        f.write(f'\n# Session end: {status}\n')
    write_sidecar(sidecar_path, args=args, messages=messages,
                  written_at_turn=turn,
                  status='Done' if is_done(last) else 'max_turns_reached',
                  starter_label=starter_label, adversary_label=adversary_label)
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
        prog='molopt_oa.py',
        description='Adversarial molecule optimization: OpenAI <-> Anthropic '
                    '(selectable starter with tools, the other critiques).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # OpenAI starts (tool-calling proposer), Anthropic critiques (default)
  python3 molopt_oa.py --protein HMGCR
  python3 molopt_oa.py --start openai --openai-model gpt-5.2

  # Anthropic starts, OpenAI critiques
  python3 molopt_oa.py --start anthropic --anthropic-model claude-haiku-4-5-20251001

  # Gemini starts (tool-calling proposer), OpenAI critiques (defaults)
  python3 molopt_oa.py --start gemini
  python3 molopt_oa.py --start gemini --adversary anthropic   # Gemini vs Anthropic

  # Quick chemistry-stack check with no LLM keys
  python3 molopt_oa.py --self-test

Keys are read from the environment (OPENAI_API_KEY, ANTHROPIC_API_KEY,
GEMINI_API_KEY) or a .env file in the working directory; CLI flags override.
Only the two providers actually used (starter + adversary) need their keys.
Source ~/.zshrc first if the keys live there.
""",
    )
    p.add_argument('--protein', default='HMGCR',
                   help='Docking target, a dockstring target name (default: HMGCR). '
                        'Docking works for any of dockstring\'s 58 targets; residue-contact '
                        'analysis (dock_and_get_interacting_residues) needs a prepared '
                        'receptor PDB on disk: HMGCR, ADRB1, ADRB2, MAOB, DRD2.')
    p.add_argument('--start', choices=['openai', 'anthropic', 'gemini'], default='openai',
                   help='Which provider leads as the tool-calling proposer; the other '
                        'becomes the critique-only adversary (default: openai).')
    p.add_argument('--adversary', choices=['openai', 'anthropic', 'gemini'], default=None,
                   help='Which provider critiques. If omitted, defaults to the natural '
                        'counterpart: openai<->anthropic pair up; gemini defaults to an '
                        'openai adversary (matches the Ollama study). May explicitly match '
                        '--start for a same-provider self-critique run (e.g. openai vs '
                        'openai) -- only the auto-derived default is required to differ.')

    p.add_argument('--openai-model', default=OPENAI_DEFAULT_MODEL,
                   help=f'OpenAI model name (default: {OPENAI_DEFAULT_MODEL}). Used as the '
                        f'starter when --start openai, else as the adversary.')
    p.add_argument('--anthropic-model', default=ANTHROPIC_DEFAULT_MODEL,
                   help=f'Anthropic model name (default: {ANTHROPIC_DEFAULT_MODEL}). Used as '
                        f'the starter when --start anthropic, else as the adversary.')
    p.add_argument('--gemini-model', default=GEMINI_DEFAULT_MODEL,
                   help=f'Gemini model name (default: {GEMINI_DEFAULT_MODEL}). Used as the '
                        f'starter when --start gemini, else as the adversary.')
    p.add_argument('--openai-key', default=None, help='OpenAI API key (or env OPENAI_API_KEY).')
    p.add_argument('--anthropic-key', default=None, help='Anthropic API key (or env ANTHROPIC_API_KEY).')
    p.add_argument('--gemini-key', default=None, help='Gemini API key (or env GEMINI_API_KEY).')
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
                        'starter turn; continues into the adversary refinement loop from '
                        'the last assistant message already in the sidecar.')
    p.add_argument('--results-dir', default=os.path.join(_HERE, 'results'),
                   help='Where to write the timestamped results .md (default: ./results).')
    p.add_argument('--max-turns', type=int, default=20, help='Safety cap on adversary<->starter turns (default: 20).')
    p.add_argument('--max-tool-calls', type=int, default=12,
                   help='Max tool-calling rounds per starter turn (default: 12). '
                        'Prevents a stuck model from looping forever; once hit, the '
                        'model is forced to emit a text response.')
    p.add_argument('--max-retries', type=int, default=3, help='Retries on API connection errors (default: 3).')
    p.add_argument('--api-timeout', type=float, default=_DEFAULT_API_TIMEOUT,
                   help=f'Per-call API timeout in seconds (default: {_DEFAULT_API_TIMEOUT:.0f}). '
                        f'Permanent errors (bad auth / out-of-credits) fail fast regardless; '
                        f'this bounds a genuinely hung connection.')
    p.add_argument('--quiet', action='store_true', help='Suppress thinking/content/tool prints.')
    p.add_argument('--rdkit-verbose', action='store_true',
                   help='Re-enable RDKit stderr (SMILES Parse Error) logs. They are silenced '
                        'by default since invalid substituents are now surfaced to the model '
                        'as "invalid SMILES, skipped" entries in the tool results.')
    p.add_argument('--trace', action='store_true',
                   help='Write the full tool-call trace into the results .md '
                        '(default: only terse section headers + text go to the md; the '
                        'trace otherwise only goes to stdout).')

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

    if args.self_test:
        self_test()
        return 0

    if args.resume:
        if not os.path.isfile(args.resume):
            raise SystemExit(f"Resume sidecar not found: {args.resume}")
    elif not os.path.isfile(args.context_file):
        raise SystemExit(f"Context file not found: {args.context_file}")

    run_session(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())