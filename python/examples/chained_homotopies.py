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

The family here is deliberately tiny but NOT all sunshine: its first solve has four
paths of which two are TRUNCATED (the endgame's security check cuts them off as they
flee toward infinity -- a verdict, not a failure), so the picture shows lineages that
die at the first column and lineages that survive the whole chain.  Rerunning the
script is (nearly) instant: every solve is ensure-answered, so recorded paths recall
instead of recomputing.
"""

import argparse

import bertini as pb
from bertini import Variable, VariableGroup, System
from bertini.nag_algorithm import blend_homotopy


X, Y = Variable('x'), Variable('y')


def family_member(a):
    """The member {x^2 - a^2, x*y - 1}: two finite roots (+-a, +-1/a) -- but a total
    degree of 4, so a fresh solve tracks four paths and two of them have nowhere
    finite to go."""
    sys = System()
    sys.add_variable_group(VariableGroup([X, Y]))
    sys.add_function(X**2 - a**2)
    sys.add_function(X * Y - 1)
    return sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', metavar='STEM', default=None,
                        help='draw the left-to-right chain progression to STEM.png and '
                             'STEM.svg (any suffix on STEM is replaced)')
    args = parser.parse_args()

    # ---- build the chain -----------------------------------------------------------
    members = [family_member(a) for a in (1, 2, 3, 4)]

    # the root: an ordinary solve (total-degree start; provenance bottoms out at
    # canonical start labels).  Four paths -- and only two finite roots, so two paths
    # get TRUNCATED near infinity: recorded as `diverged`, an answer in itself.
    results = [pb.solve(members[0], seed=42)]

    # each further member: continue the previous solutions through a blend homotopy.
    # Only the FINITE solutions are carried forward -- the diverged paths' lineages
    # simply end, which is exactly what the picture shows.
    for previous, target in zip(members, members[1:]):
        homotopy = blend_homotopy(target, previous)
        results.append(pb.solve(target, homotopy=homotopy, start=results[-1], seed=42))

    # margin notes and deliverables travel with the records
    final = results[-1]
    pb.annotate(final.solutions[0], 'note', 'the branch through x = +1')
    pb.save('the chained family', final,
            description='solutions of the last member, chained from the first')

    # ---- read the records back: the navigation tools --------------------------------
    print('the runs, one row each:\n')
    print(pb.runs()[['run', 'when', 'num_paths', 'seed']].to_string(index=False))

    print('\nthe tracked paths (statuses and where each one started):\n')
    tracks = pb.tracks()
    print(tracks[['run', 'index', 'status', 'outcome', 'start_kind']].to_string(index=False))
    # the first solve's two extra total-degree paths were truncated near infinity;
    # everything that had somewhere to go, went there
    assert (tracks['status'] == 'diverged').sum() == 2
    assert (tracks['status'] == 'success').sum() == len(tracks) - 2

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
        from pathlib import Path
        ax = pb.plot_chain()
        stem = Path(args.plot).with_suffix('')      # tutorial images ship as png AND svg
        for suffix in ('.png', '.svg'):
            ax.figure.savefig(stem.with_suffix(suffix), dpi=110, bbox_inches='tight')
        print('\nchain progression drawn to %s.png / .svg' % stem)

    print('\nrecords at:', pb.records_dir(), '-- results.json has the story.')


if __name__ == '__main__':
    main()
