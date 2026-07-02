# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""The kill-and-rerun demo: the arc's acceptance test, live.

    PYTHONPATH=python python prototypes/ledger_v0/demo_resume.py [ledger_dir]

Solves a degree-24 system (2*3*4 start points), dies a simulated walltime death partway
through, then re-invokes the SAME ensure_solved call, which resumes from the journal and
finishes.  Poke the ledger afterward with nothing but standard tools:

    jq .kind    <ledger>/journals/*.jsonl | sort | uniq -c
    cat <ledger>/objects/*/*            # the target system, readable classic input
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable

from ledger import Ledger
from memo_solve import ensure_solved, SimulatedCrash


def build_target():
    x, y, z = Variable("x"), Variable("y"), Variable("z")
    s = System()
    s.add_variable_group(VariableGroup([x, y, z]))
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x**3 - z)
    s.add_function(y**4 + z**2 - 2)
    return s


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "demo_ledger"
    ledger = Ledger(root)
    print("ledger at:", Path(root).resolve())

    print("\n--- run 1: dies after 9 of 24 paths (simulated walltime kill) ---")
    try:
        ensure_solved(build_target(), ledger, crash_after=9)
    except SimulatedCrash as crash:
        print("CRASH:", crash)

    print("\n--- run 2: the same call again (this is 'resume': there is no resume) ---")
    result = ensure_solved(build_target(), ledger)
    print("reused from ledger: %d paths" % result.num_reused)
    print("computed now:       %d paths" % result.num_computed)
    print("total finite endpoints recorded: %d" % len(result.solutions))

    print("\n--- run 3: rerun of the completed ask is a no-op ---")
    result = ensure_solved(build_target(), ledger)
    print("reused: %d   computed: %d" % (result.num_reused, result.num_computed))


    print("\n--- what is on disk ---")
    print(" ", ledger.describe())


if __name__ == "__main__":
    main()
