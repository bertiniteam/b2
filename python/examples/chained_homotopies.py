#!/usr/bin/env python3
"""Chained homotopies with full provenance: the records follow every hop.

A *chain* is a sequence of solves where each solve's start points are the previous
solve's solutions -- the workhorse pattern of parameter continuation.  With the
structured output directory attached (which ``bertini.solve`` does by default), the
chain is *recorded*: every tracked path remembers which point it started from, so any
final solution can be walked back through every intermediate run to the very first
start point.

This script builds a small chain, then reads its own records back with the navigation
tools -- the runs and tracks as pandas DataFrames, the provenance as a networkx graph,
and the whole chain drawn left to right::

    python chained_homotopies.py                       # writes ./bertini_output
    python chained_homotopies.py --plot chain.png      # also draw the progression

The family here is deliberately tiny -- circles of growing radius intersected with the
line x = y -- so the chain itself is the star.  Rerunning the script is (nearly)
instant: every solve is ensure-answered, so recorded paths hydrate instead of
recomputing.
"""

import argparse

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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', metavar='PATH', default=None,
                        help='draw the left-to-right chain progression to this file')
    args = parser.parse_args()

    # ---- build the chain -----------------------------------------------------------
    radii_squared = [1, 4, 9, 16]
    members = [circle_line(r2) for r2 in radii_squared]

    # the root: an ordinary solve (total-degree start; provenance bottoms out at
    # canonical start labels)
    results = [pb.solve(members[0], seed=42)]

    # each further member: continue the previous solutions through a blend homotopy
    for previous, target in zip(members, members[1:]):
        homotopy = blend_homotopy(target, previous)
        results.append(pb.solve(target, homotopy=homotopy, start=results[-1], seed=42))

    # margin notes and deliverables travel with the records
    final = results[-1]
    pb.annotate(final.solutions[0], 'note', 'the positive branch at r=4')
    pb.save('the chained family', final,
            description='solutions of the last member, chained from the first')

    # ---- read the records back: the navigation tools --------------------------------
    print('the runs, one row each:\n')
    print(pb.runs()[['run', 'when', 'num_paths', 'seed']].to_string(index=False))

    print('\nthe tracked paths (statuses and where each one started):\n')
    tracks = pb.tracks()
    print(tracks[['run', 'index', 'status', 'start_kind']].to_string(index=False))
    assert (tracks['status'] == 'success').all()

    # the provenance graph: every path an edge from its start to its endpoint
    graph = pb.provenance_graph()
    print('\nprovenance graph: %d points, %d edges'
          % (graph.number_of_nodes(), graph.number_of_edges()))

    # the walk itself: one final solution, back to the very beginning
    trail = pb.provenance(final.solutions[0])
    print('\nprovenance of one final solution (newest hop first):')
    for hop in trail:
        print('   ', hop)
    assert trail[-1]['kind'] in ('start_label', 'given_ref')

    # ---- the picture: paths flowing left to right through the chain -----------------
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        ax = pb.plot_chain()
        ax.figure.savefig(args.plot, dpi=110, bbox_inches='tight')
        print('\nchain progression drawn to', args.plot)

    print('\nrecords at:', pb.records_dir(), '-- results.json has the story.')


if __name__ == '__main__':
    main()
