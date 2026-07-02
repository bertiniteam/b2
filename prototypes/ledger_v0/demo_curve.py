# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""The long-chain demo: a numerical real cellular decomposition SKETCH for a plane
curve, done entirely as memoized, recorded homotopy runs.

    PYTHONPATH=python python prototypes/ledger_v0/demo_curve.py [dir]

Curve: the ellipse x^2 + 4 y^2 = 4, projected to the x-axis.

  stage 1  critical points of the projection: solve { f, df/dy }        (base run)
  stage 2  witness slice of the family S_c = { f, x - c } at generic c  (base run)
  stage 3  midpoint slice of the bounded interval                        (chain, depth 2)
  stage 4  SAMPLE each edge at several x values                          (chains, depth 3)
           -- each sample point annotated with its projection value.

Run 1 dies a simulated walltime death during sampling; run 2 (the same script) resumes;
then one sample point's provenance is walked back four links to its total-degree start
label.  This is bertini_real's curve-case shape in miniature: the chain-of-homotopies
demo that NID would otherwise provide.  (Edge membership here is inferred from the sign
of y -- a sketch, not the real connect-the-dots machinery.)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from bertini import System, Variable, VariableGroup
from bertini.symbolics import Rational

from ledger import Ledger
from memo_solve import (solve, continue_from, provenance_chain,
                        annotate, annotations_for, declare_result, SimulatedCrash)


def curve_f(x, y):
    return x**2 + 4 * y**2 - 4


def critical_system():
    """{ f, df/dy }: the projection-critical points of the ellipse."""
    x, y = Variable("x"), Variable("y")
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(curve_f(x, y))
    s.add_function(8 * y)                          # df/dy
    return s


def slice_at(c_num, c_den=1):
    """S_c = { f, x - c } with c exact rational (identity-stable across reruns)."""
    x, y = Variable("x"), Variable("y")
    s = System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(curve_f(x, y))
    s.add_function(x - Rational("%d/%d" % (c_num, c_den)))
    return s


# sampling plan: (numerator, denominator) x-values inside the interval (-2, 2)
SAMPLE_XS = [(-3, 2), (-1, 2), (1, 2), (3, 2)]


def decompose(ledger, crash_after=None):
    """The whole decomposition as one rerunnable script (ensure-answered throughout)."""
    print("  stage 1: critical points of the projection")
    crit = solve(critical_system(), ledger)
    crit_xs = sorted(complex(sol[0]).real for sol in crit.solutions.values()
                     if abs(complex(sol[0]).imag) < 1e-8)
    print("    critical x: %s   (reused %d, computed %d)"
          % (crit_xs, crit.num_reused, crit.num_computed))

    print("  stage 2: witness slice at generic c = 1/3")
    witness = solve(slice_at(1, 3), ledger)
    print("    %d points on the slice  (reused %d, computed %d)"
          % (len(witness.solutions), witness.num_reused, witness.num_computed))

    print("  stage 3: midpoint slice of the bounded interval (c = 0)")
    midpoint = continue_from(slice_at(0), slice_at(1, 3), ledger)
    edges = {i: ("top" if complex(sol[1]).real > 0 else "bottom")
             for i, sol in midpoint.solutions.items()}
    print("    edges: %s  (reused %d, computed %d)"
          % (edges, midpoint.num_reused, midpoint.num_computed))

    print("  stage 4: sampling the edges")
    sample_runs = []
    for num, den in SAMPLE_XS:
        result = continue_from(slice_at(num, den), slice_at(0), ledger,
                                  crash_after=crash_after)
        xval = num / den
        for i, sol in result.solutions.items():
            if not annotations_for(ledger, result.run_id, i):
                annotate(ledger, result.run_id, i,
                         projection=xval, edge=edges.get(i, "?"))
        print("    x = %+.1f: %d points  (reused %d, computed %d)"
              % (xval, len(result.solutions), result.num_reused, result.num_computed))
        sample_runs.append(result)

    # the deliverable: only the sample points are RESULTS; everything else was scaffolding
    declare_result(ledger, "edge samples of the ellipse x^2+4y^2=4",
                   [(r.run_id, i) for r in sample_runs for i in sorted(r.solutions)],
                   description="x-projection cellular decomposition sketch; "
                               "2 edges sampled at 4 x-values each")
    return sample_runs


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "demo_curve_ledger"
    ledger = Ledger(root)
    print("records at:", Path(root).resolve())

    print("\n--- run 1: dies during edge sampling ---")
    try:
        decompose(ledger, crash_after=1)
    except SimulatedCrash as crash:
        print("  CRASH:", crash)

    print("\n--- run 2: the same script; everything done is skipped ---")
    samples = decompose(ledger)

    print("\n--- a sample point, annotated: ---")
    last = samples[-1]
    print(" ", annotations_for(ledger, last.run_id, 0))

    print("\n--- its provenance, back to the beginning (4 links) ---")
    for link in provenance_chain(ledger, last.run_id, 0):
        print(" ", link)


    print("\n--- what is on disk ---")
    print(" ", ledger.describe())


if __name__ == "__main__":
    main()
