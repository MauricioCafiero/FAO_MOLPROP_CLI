#!/bin/zsh
# run_orca_validation.sh - Sequentially run the prepared ORCA singlepoints in
# results/batches/hl_batches/orca_validation/agentic_5x4_top/ (one molecule at
# a time, 4 procs / 4 GB each), skipping any molecule whose .out already ends
# with a final single-point energy (resumable).
#
# Launch detached, NOT via a harness-tracked background task (this mac reaps
# those after ~30-50 min):
#   ~/miniforge3/bin/python -c "import subprocess; subprocess.Popen(['zsh','.../run_orca_validation.sh'], stdout=open('.../orca_validation/run.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, cwd='.../orca_validation')"
# Verify with: ps -ax -o pid,ppid,sess,command | grep run_orca_validation
# Monitor by tailing run.log; do NOT sleep-poll.

caffeinate -w $$ &   # hold off idle sleep; dies with this script

ORCA=/Users/cafierom/opt/orca-6.1.1/orca
export PATH="/Users/cafierom/opt/openmpi-4.1.8/bin:$PATH"
export OMP_NUM_THREADS=4

ROOT="$(cd "$(dirname "$0")" && pwd)/.."   # FAO_MOLPROP_CLI root (script lives in code_hl/)
TOP="$ROOT/results/batches/hl_batches/orca_validation/agentic_5x4_top"
cd "$TOP" || exit 1

for d in $TOP/*/; do
  slug=$(basename "$d")
  [ -d "$d" ] || continue
  inp=$(ls "$d"*.inp 2>/dev/null | head -1)
  [ -n "$inp" ] || continue
  out="${inp%.inp}.out"
  if grep -q "FINAL SINGLE POINT ENERGY" "$out" 2>/dev/null; then
    echo "[skip] $slug already done"
    continue
  fi
  echo "[run ] $slug  $(date '+%H:%M:%S')"
  (cd "$d" && "$ORCA" "$(basename "$inp")" > "$(basename "$out")" 2>&1)
  if grep -q "FINAL SINGLE POINT ENERGY" "$out"; then
    echo "[done] $slug  $(date '+%H:%M:%S')"
  else
    echo "[FAIL] $slug  $(date '+%H:%M:%S')  (see $out)"
  fi
done
echo "all molecules processed  $(date)"