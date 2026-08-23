#!/usr/bin/env python3
"""One-shot smoke test for GeminiActor.chat_turn (thought_signature round-trip).

Runs a single proposer turn that forces a calculate_SAS_and_NP tool call, then
verifies the tool loop completed without the 400 thought_signature error and
that signatures were preserved into the stored message dicts. CPU-light: RDKit
only, no Vina docking. Prints a final PASS/FAIL line.
"""
import os, sys
sys.argv = ['x']
import molopt_oa as m

key = os.environ.get('GEMINI_API_KEY')
assert key, "GEMINI_API_KEY missing -- source ~/.zshrc"

actor = m.GeminiActor(model='gemini-3-flash-preview', api_key=key)
sys_msg = m.build_system_message()
messages = [
    {'role': 'system', 'content': sys_msg},
    {'role': 'user', 'parts': [{'text': (
        'Here is my starting molecule list:\n'
        '- c1ccccc1 (benzene), score -5.0\n'
        '- c1ccc(O)cc1 (phenol), score -5.5')}]},
]
prompt = ("Propose ONE improved analogue and BEFORE you finalize it, call "
          "calculate_SAS_and_NP on your candidate to check its synthetic "
          "accessibility. Then give your final proposal with the SAS/NP scores.")

try:
    msgs, last, trace = actor.chat_turn(
        messages, prompt, max_retries=2, max_tool_calls=4, verbose=True)
except Exception as err:
    print(f"\n=== SMOKE TEST FAILED (exception) ===\n{err!r}")
    sys.exit(1)

print("\n=== CHAT_TURN RETURNED (no 400 thought_signature error) ===")
print("last_assistant_text:", (last or '')[:400])
print("n messages:", len(msgs))
n_sig = 0
for i, mm in enumerate(msgs):
    if mm.get('role') in ('model', 'user'):
        for p in mm.get('parts', []):
            if 'thought_signature' in p and p['thought_signature']:
                n_sig += 1
                print(f"  msg[{i}] {mm['role']}: part carries thought_signature "
                      f"(b64 len {len(p['thought_signature'])}) keys={list(p.keys())}")
print(f"parts with preserved thought_signature: {n_sig}")
# A successful tool round-trip must have at least the model function_call sig AND
# the echoed user function_response sig (>=2).
ok = n_sig >= 2 and bool(last)
print("=== SMOKE TEST PASSED ===" if ok else "=== SMOKE TEST FAILED (no sigs / no text) ===")
sys.exit(0 if ok else 1)