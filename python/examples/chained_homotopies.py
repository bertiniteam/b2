#!/usr/bin/env python3
"""Chained homotopies with full provenance: the records follow every hop.

A *chain* is a sequence of solves where each solve's start points are the previous
solve's solutions -- the workhorse pattern of parameter continuation.  With the
structured output directory attached (which ``bertini.solve`` does by default), the
chain is *recorded*: every tracked path remembers which point it started from, so any
final solution can be walked back through every intermediate run to the very first
start point.  That walk is what this example demonstrates::

    python chained_homotopies.py            # writes ./bertini_output (or set
                                            # BERTINI_RECORDS_DIR)

The family here is deliberately tiny -- circles of growing radius intersected with the
line x = y -- so the chain itself is the star:

  1. solve the first family member from scratch (a total-degree start);
  2. for each next member, build the gamma-trick blend homotopy from the previous
     member and continue the previous solutions to the new ones;
  3. annotate a point, declare the deliverables, and walk the provenance.

Rerunning this script is (nearly) instant: every solve is ensure-answered, so the
recorded paths hydrate instead of recomputing.
"""

import bertini as pb
from bertini import Variable, VariableGroup, System
from bertini.nag_algorithm import blend_homotopy


def circle_line(radius_squared):
    """The circle x^2 + y^2 = r^2 intersected with the line x = y."""
    x, y = Variable('x'), Variable('y')
    sys = System()
    sys.add_variable_group(VariableGroup([x, y]))
    sys.add_function(x**2 + y**2 - radius_squared)
    sys.add_function(x - y)
    return sys


def main():
    radii_squared = [1, 4, 9, 16]
    members = [circle_line(r2) for r2 in radii_squared]

    # the root of the chain: an ordinary solve (total-degree start; provenance
    # bottoms out at canonical start labels)
    results = [pb.solve(members[0], seed=42)]

    # each further member: continue the previous solutions through a blend homotopy
    for previous, target in zip(members, members[1:]):
        homotopy = blend_homotopy(target, previous)
        results.append(pb.solve(target, homotopy=homotopy, start=results[-1], seed=42))

    for r2, result in zip(radii_squared, results):
        print("r^2 = %-3d run %s: %d solutions, %d hydrated"
              % (r2, result.run_id, len(result), result.num_hydrated))

    # margin notes travel with the records
    final = results[-1]
    pb.annotate(final.solutions[0], 'note', 'the positive branch at r=4')
    pb.save('the chained family', final,
            description='solutions of the last member, chained from the first')

    # the dream, demonstrated: walk a final point back to the very beginning
    trail = pb.provenance(final.solutions[0])
    print("\nprovenance of one final solution (newest hop first):")
    for hop in trail:
        print("   ", hop)
    assert trail[-1]['kind'] in ('start_label', 'given_ref')

    print("\nrecords at:", pb.records_dir(), "-- results.json has the story.")


if __name__ == '__main__':
    main()
