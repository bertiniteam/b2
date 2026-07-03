Chained homotopies: provenance all the way back
###############################################

A *chain* is a sequence of solves where each solve's start points are the previous
solve's solutions — the workhorse pattern of parameter continuation.  Because every
solve records where each of its paths started (see :doc:`../your_records/index`), a
chain is more than its final answers: it is a **directed graph of points**, and any
final solution can be walked back through every intermediate run to the very first
start point.

This tutorial builds a small chain — circles of growing radius intersected with the
line :math:`x = y` — and then reads its own records back with the navigation tools.
The complete script is :download:`chained_homotopies.py
<../../../../examples/chained_homotopies.py>`:

.. literalinclude:: ../../../../examples/chained_homotopies.py
   :language: python
   :caption: examples/chained_homotopies.py

Chaining is two keyword arguments
=================================

The root of a chain is an ordinary solve.  Every further link passes the homotopy you
built and where its paths start:

.. code-block:: python

   results.append(pb.solve(target, homotopy=blend_homotopy(target, previous),
                           start=results[-1], seed=42))

``start=`` accepts a prior :class:`~bertini.records.SolveResult` (or its solutions) —
those points carry provenance, so the new run's records link back to them with
``point_ref`` references.  Raw arrays work too: they are archived as a **given**, and
provenance bottoms out honestly at data you supplied.  To chain from a run recorded in
*another session* — or by the command-line ``bertini2`` — read its endpoints cold with
:func:`bertini.solutions_of`; they come back as points with provenance, ready to chain.

Reading the chain back
======================

The navigation tools read the plain records — any directory, any producer, no solver
objects:

* :func:`bertini.runs` — one :class:`pandas.DataFrame` row per run: when, how many
  paths, which seed, which software wrote it.
* :func:`bertini.tracks` — one row per tracked path: its verdict (``success`` /
  ``diverged`` / ``failed``), and where it started.  Endpoint coordinates stay out of
  the table unless you pass ``coordinates=True`` — on a million-path audit you want
  the statuses, not the bytes.
* :func:`bertini.provenance_graph` — the points as a :class:`networkx.DiGraph`, one
  edge per path, pointing from where it started to where it ended.  Time flows along
  the edges.
* :func:`bertini.provenance` — the walk for a single point: its chain of
  ``{'run', 'index'}`` hops, ending at a start label or a given.

The picture
===========

:func:`bertini.plot_chain` draws the chain **left to right** — each run a column,
paths flowing rightward from their origins (squares) through every run, colored by
verdict (green success, orange diverged, red failed):

.. image:: chain_progression.png
   :alt: paths flowing left to right through a four-run chain
   :align: center

Large solves are a design constraint, not an afterthought: a run with more than
``max_paths_drawn`` paths (default 200) is not drawn path-by-path — the figure
automatically aggregates to one node per run, with edge widths showing how many paths
flow between runs, so a million-path chain renders instead of crashing your session.
The same instinct applies to :func:`~bertini.provenance_graph`: pass ``runs=`` to
restrict a big directory to the chains you care about before building a graph out of
it.

Regenerating the figure
=======================

The image above is committed; regenerate it after changing the example or the plotting
code::

   cd python/docs/source/tutorials/chained_homotopies
   BERTINI_RECORDS_DIR=$(mktemp -d) python ../../../../examples/chained_homotopies.py \
       --plot chain_progression.png
