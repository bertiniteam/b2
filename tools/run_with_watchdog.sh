#!/usr/bin/env bash
# run_with_watchdog.sh -- run a command with a hard wall-clock timeout, killing it (and its whole
# process group) if it exceeds the limit.
#
# WHY THIS EXISTS: some Bertini 2 solves can spin or churn for a very long time (notably the
# linear-product TotalDegree start system, whose Cauchy endgame can stall near t=0 -- see
# z_notes/20260629_endgame_stepsize_reset_rootcause).  Running such a solve directly from an agent
# session (or CI) can peg a core indefinitely and starve everything else.  macOS has no `timeout(1)`,
# so use this wrapper for ANY command that might hang -- especially `... python -c 'zd.solve()'` and
# the `test_*` binaries.
#
# Usage:
#   tools/run_with_watchdog.sh <seconds> <command> [args...]
# Examples:
#   tools/run_with_watchdog.sh 30 ./build/core/test_nag_algorithms
#   OMP_NUM_THREADS=1 PYTHONPATH=python tools/run_with_watchdog.sh 60 python3 my_solve.py
#
# Exit code: the command's exit code if it finished in time; 124 if it was killed by the watchdog
# (the conventional timeout(1) code).
set -u

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <seconds> <command> [args...]" >&2
  exit 2
fi

timeout_s="$1"; shift

# New process group so we can kill the command AND any children (worker threads/processes) it spawns.
set -m
"$@" &
cmd_pid=$!

elapsed=0
while kill -0 "$cmd_pid" 2>/dev/null; do
  sleep 1
  elapsed=$((elapsed + 1))
  if [ "$elapsed" -ge "$timeout_s" ]; then
    echo "[watchdog] '$*' exceeded ${timeout_s}s -- killing (likely a hang/churn)." >&2
    kill -9 -- "-${cmd_pid}" 2>/dev/null || kill -9 "$cmd_pid" 2>/dev/null
    wait "$cmd_pid" 2>/dev/null
    exit 124
  fi
done

wait "$cmd_pid"
exit $?
