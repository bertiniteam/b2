# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""The chain demo: a parameter sweep as a fan of two-link provenance chains.

    PYTHONPATH=python python prototypes/ledger_v0/demo_chain.py [dir]

Solves one generic member of the family { x^2 + y^2 = a, x = y } from scratch (the
expensive ancestor), then sweeps four targets by parameter continuation -- dying a
simulated walltime death mid-sweep, resuming with the same calls, and finally walking
one endpoint's provenance all the way back to its total-degree start label.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from bertini import System, VariableGroup
from bertini.function_tree.symbol import Variable

from ledger import Ledger
from memo_solve import solve, continue_from, provenance_chain, SimulatedCrash


def circle_family(a):
    x, y = Variable("x"), Variable("y")
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - a)
    s.add_function(x - y)
    return s


SWEEP = [2, 3, 5, 7]


def sweep(ledger, crash_after_total=None):
    """The user's script: solve the generic, continue to each target.  Identical on
    first run, after a crash, and after completion -- ensure-answered all the way."""
    generic = circle_family(13)
    base = solve(generic, ledger)
    print("  generic a=13: reused %d, computed %d" % (base.num_reused, base.num_computed))

    budget = crash_after_total
    for a in SWEEP:
        result = continue_from(circle_family(a), generic, ledger,
                                  crash_after=budget)
        print("  target a=%d: reused %d, computed %d" % (a, result.num_reused, result.num_computed))
        if budget is not None:
            budget -= result.num_computed
            if budget <= 0:
                budget = None
    return base


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "demo_chain_ledger"
    ledger = Ledger(root)
    print("records at:", Path(root).resolve())

    print("\n--- run 1: the sweep dies mid-flight (after 3 continued paths total) ---")
    try:
        sweep(ledger, crash_after_total=3)
    except SimulatedCrash as crash:
        print("  CRASH:", crash)

    print("\n--- run 2: the same script again; finished work is skipped ---")
    base = sweep(ledger)

    print("\n--- run 3: no-op ---")
    sweep(ledger)

    print("\n--- provenance: one endpoint of the a=7 solve, back to the beginning ---")
    last = continue_from(circle_family(7), circle_family(13), ledger)
    for link in provenance_chain(ledger, last.run_id, 0):
        print(" ", link)


    print("\n--- what is on disk ---")
    print(" ", ledger.describe())


if __name__ == "__main__":
    main()
