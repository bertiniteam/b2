Your results, your records
===========================

Every solve writes a **structured output directory**: a durable, plain-text record of
what was computed — which system, which settings, which seed, and every path's endpoint
with its provenance.  You are free to delete it (the only consequence is recomputing),
free to ``grep`` it, and free to read it in ten years with no bertini installed: the
directory carries its own ``README.txt`` explaining the format.

How to think about it: the records are ordinary program *output*, like the ``.log``
file LaTeX writes — always produced, never precious, occasionally exactly what you
need.  They are **not** a database you administer.  There is nothing to configure, no
schema to migrate, no server; just files you can read.

What it gains you:

* **Restartability.**  A solve is *ensure-answered*: it consults the records first,
  takes any paths already answered for the identical ask (system + settings + seed),
  and computes only the rest.  Kill a 40-hour run at hour 39 and rerunning the same
  command finishes the last hour.  This works for the CLI and Python alike, with no
  flags — resuming *is* rerunning.
* **Reproducibility.**  ``seed=42`` names the exact homotopy — same gamma, start
  points, and patch, on every machine, forever.  Even a run where you *omitted* the
  seed records the one it drew, so accidental results remain reproducible.
* **Provenance.**  Every path remembers where it started; chained solves link runs
  into a walkable graph (see :doc:`../chained_homotopies/index`), so any final point
  answers "which chain of homotopies and start points produced you?" — all the way to
  the beginning.
* **An audit trail.**  Every tracked path is recorded whatever its outcome —
  ``success``, ``diverged`` (a path that went to infinity: an answer, not a failure),
  or ``failed`` (the tracker gave up) — so "what happened to my 50 missing paths?" is
  a one-line pandas query, not a mystery.

Turning it on and off
---------------------

It is on by default, everywhere.  Off is one line:

.. testcode::

   import bertini as pb
   pb.recording(False)     # solves run bare: nothing written, nothing consulted
   pb.recording(True)      # back on
   print(pb.recording())   # ask the current state

.. testoutput::

   True

For the command line, the switch is the environment: ``BERTINI_RECORDS_DIR=""``
(empty) runs ``bertini2`` with no records, and any non-empty value relocates them
(``BERTINI_RECORDS_DIR=~/project/records bertini2 input``).  With recording off there
is nothing to resume from — the trade is yours to make.

The solver classes (:class:`~bertini.ZeroDimSolver`, :class:`~bertini.HomotopySolver`)
record *ambiently*: they join in once the ambient directory has been named, either by
``bertini.records_dir("my_records")`` in your script (one line — see
:doc:`../parameter_homotopy/index` for it in action) or by the ``BERTINI_RECORDS_DIR``
environment variable.  Until then they run bare, exactly as they always have.

The three verbs
---------------

.. testcode::

   import bertini as pb

   x, y = pb.Variable('x'), pb.Variable('y')
   s = pb.System()
   s.add_variable_group(pb.VariableGroup([x, y]))
   s.add_function(x**2 + y**2 - 1)
   s.add_function(x - y)

   import tempfile, os
   where = os.path.join(tempfile.mkdtemp(), 'records')   # a scratch spot for this tutorial;
                                                         # normally you write NOTHING and get ./bertini_output

   sols = pb.solve(s, seed=42, directory=where)
   print(len(sols), 'solutions, recalled', sols.num_recalled)

.. testoutput::

   2 solutions, recalled 0

Solve it *again* — same system, same seed — and nothing is recomputed: the answers come
back from the records.

.. testcode::

   again = pb.solve(s, seed=42, directory=where)
   print('recalled', again.num_recalled, 'of', len(again))

.. testoutput::

   recalled 2 of 2

Each solution is the numpy-like point you expect, and it *remembers where it came from*:

.. testcode::

   pt = sols[0]
   print(sorted(pt.provenance))

.. testoutput::

   ['index', 'run']

``save`` anything under a name; ``load`` it back — in this session or any later one:

.. testcode::

   pb.save('my favorites', sols, directory=where)
   pb.save('notes', {'count': 2}, directory=where)
   print(pb.load('notes', directory=where)['value'])

.. testoutput::

   {'count': 2}

Seeds are identities
--------------------

``seed=42`` means the *same homotopy* — the same gamma, start points, and patch —
forever, on every machine.  A seed you liked is a run you can reproduce, share, and
resume.  Omit the seed and a fresh one is drawn (and recorded, so even accidental runs
are reproducible afterward).

Chains: provenance all the way back
-----------------------------------

A *chained* solve continues a previous result's solutions through a homotopy you
build — and the records follow every hop.  Chain with ``homotopy=`` and ``start=``:

.. testcode::

   from bertini.nag_algorithm import blend_homotopy

   bigger = pb.System()
   bigger.add_variable_group(pb.VariableGroup([x, y]))
   bigger.add_function(x**2 + y**2 - 4)     # the same family, radius 2
   bigger.add_function(x - y)

   H = blend_homotopy(bigger, s)            # from the circle we already solved
   sols2 = pb.solve(bigger, homotopy=H, start=sols, seed=42, directory=where)
   print(len(sols2))

.. testoutput::

   2

Every new endpoint remembers which point it came from.  ``provenance`` walks the
chain back to the very beginning — through as many runs as it takes, including runs
written by the command-line ``bertini2``:

.. testcode::

   trail = pb.provenance(sols2.solutions[0], directory=where)
   print(trail[0]['run'] == sols2.run_id)
   print(trail[-1]['kind'])

.. testoutput::

   True
   start_label

``solutions_of`` reads any run's endpoints straight from the records — no solver
object, any session, any machine — as points that chain directly:

.. testcode::

   cold = pb.solutions_of(sols.run_id, directory=where)
   print(len(cold), cold[0].provenance['run'] == sols.run_id)

.. testoutput::

   2 True

Start points that are *your* data (arrays, not a prior result) are archived as a
**given**: provenance bottoms out honestly at what you supplied.  And margin notes
travel with the points:

.. testcode::

   pb.annotate(sols2.solutions[0], 'note', 'the positive branch', directory=where)

A complete runnable chain lives in ``python/examples/chained_homotopies.py``.

What is in the directory
------------------------

Nothing here needs bertini to read::

   README.txt      what this is, and the record format — self-contained
   results.json    what you saved, pretty-printed, plus references to the system
                   and configs that produced it: one json.load away
   INDEX.txt       one line per run: when, what, how many paths
   history/        every record, one JSON object per line (grep / jq / pandas)
   definitions/    what the records refer to, grouped by kind and named by their
                   digest -- systems/ (JSON embedding the exact canonical encoding,
                   whose sha256 IS the system's identity digest), configs/ (JSON),
                   givens/ (your data, byte-exact)

The command line gets the same treatment: running ``bertini2`` on an input file writes
``bertini_output`` beside your familiar Bertini 1 files (``main_data``,
``finite_solutions``, ...), and a killed run *finishes* when you simply run the same
command again.  Point the records somewhere else — a project directory on a cluster,
never scratch — with the ``BERTINI_RECORDS_DIR`` environment variable (under MPI, pass
``mpirun -x BERTINI_RECORDS_DIR``; only the manager rank writes).

Navigating what you have
------------------------

Four tools read the plain records back — any directory, any producer, no solver
objects: :func:`bertini.runs` and :func:`bertini.tracks` (pandas DataFrames of the
runs and the tracked paths), :func:`bertini.provenance_graph` (a networkx graph of
the points), and :func:`bertini.plot_chain` (the chain drawn left to right).  They
are the subject of :doc:`../chained_homotopies/index`.

Power users: the solver objects underneath expose the same machinery —
``solver.record_to(path)``, ``solver.num_paths_recalled()``, ``solver.records_run_id()``
— and the full record schema is documented in ``docs/records/ledgerrec-1.md``.
