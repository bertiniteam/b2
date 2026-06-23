🚀 Solving at scale: many cores, many machines
*************************************************

Homotopy continuation is *pleasantly parallel*.  Every path is independent -- once the start
points are fixed, tracking solution #3 needs nothing from solution #57.  So the wall-clock cost of
a solve is, in principle, the cost of the slowest single path, no matter how many paths there are,
*if* you have enough workers to go around.

Bertini turns that "in principle" into a one-line change.  This tutorial shows the two knobs:

* **MPI ranks** -- a *manager-worker* pool, possibly spread across many machines.  Pass a
  communicator to :func:`solve` and rank 0 hands paths out to the rest.
* **threads per rank** -- each worker can track several paths at once.  Set ``OMP_NUM_THREADS``.

The two scale *multiplicatively*: ``R`` worker ranks of ``T`` threads each track ``R × T`` paths
at a time.  We demonstrate on two driving problems -- a **cyclic-n** system (hundreds of cheap
paths) and an **eigenvalue** solve (a few dozen heavy, uneven paths) -- because they stress the
scheduler in opposite ways, and the same manager-worker pool handles both.

.. note::

   The timing tables below were **measured on an 8-core machine**.  Your numbers will differ; the
   *shape* (speedup grows with workers, then flattens once the single slowest path dominates the
   wall) is the point.  Run ``tools/update_scaling_timings.py`` to repopulate every table for your
   own hardware -- it pins a fixed seed so every layout solves the identical problem.

The manager-worker model
========================

A bertini solve already has a definite shape -- track every path *before* the endgame, run a
midpath check, then run the endgame on each path -- and the parallel solve keeps **exactly** that
shape.  There is no separate "parallel solver": you call the same :func:`solve`, only with a
communicator.

* **Rank 0 is the manager.**  It owns the bookkeeping -- the list of path indices still to do, and
  the results as they come back -- but it does not track paths itself.  It hands out work on
  demand: when a worker reports a finished path, the manager sends it the next index.  That
  *demand-driven* dispatch is what gives automatic load balancing -- a worker that draws a cheap
  path simply comes back for more sooner.
* **Ranks 1..N-1 are workers.**  Each asks the manager for a path index, tracks that path through
  the whole pre-endgame / midpath / endgame pipeline, ships the result back, and asks again, until
  the manager says "done".

Because the manager does not track, a distributed run needs **at least two ranks** (one manager +
one worker).  The serial baseline is just :func:`solve` with no communicator.  The minimal program
is::

    from mpi4py import MPI
    import bertini as pb

    comm = MPI.COMM_WORLD
    solver = pb.nag_algorithm.ZeroDim(my_system, mptype='double')

    if comm.Get_size() > 1:
        solver.solve(communicator=comm)      # rank 0 dispatches; the others track
    else:
        solver.solve()                       # serial: a lone manager would have no workers

    if pb.parallel.is_manager():             # only rank 0 has the collected solutions
        print(len(solver.solutions()), "solutions")

The :mod:`bertini.parallel` helpers -- :func:`~bertini.parallel.rank`,
:func:`~bertini.parallel.size`, :func:`~bertini.parallel.is_manager`,
:func:`~bertini.parallel.is_worker` -- let each rank know its role.  Only the manager ends up
holding the solution list, so guard your reporting and your correctness checks with
``is_manager()``.

Launch it with ``mpirun``::

    python            my_solve.py     # serial baseline (1 process)
    mpirun -n 2 python my_solve.py    # 1 manager + 1 worker  (the smallest parallel run)
    mpirun -n 9 python my_solve.py    # 1 manager + 8 workers

.. warning::

   ``mpirun`` gives stdin only to rank 0, so a heredoc (``mpirun -n 4 python << EOF``) **deadlocks**
   -- the workers wait forever for a script they never receive.  Always run a **file**.

.. note::

   To put one worker on each of your cores you want ``-n (cores + 1)`` -- the extra rank is the
   near-idle manager.  On a machine with exactly that many cores ``mpirun`` will refuse to start
   ("not enough slots"); add ``--map-by :OVERSUBSCRIBE`` to let the manager share a core (it barely
   uses one).  For example, to run 8 workers on an 8-core box::

       mpirun -n 9 --map-by :OVERSUBSCRIBE python my_solve.py

Scaling a cyclic system
=======================

The cyclic-:math:`n` system is a classic benchmark: :math:`n` polynomials whose total-degree
homotopy tracks :math:`n!` paths, most of which diverge, leaving a known number of finite
solutions.  It is the "many cheap paths" regime -- cyclic-6 is 720 paths -- so it shows off
fine-grained load balancing.

The complete, runnable script is :download:`solve_cyclic.py
<../../../examples/solve_cyclic.py>`:

.. literalinclude:: ../../../examples/solve_cyclic.py
   :language: python
   :caption: examples/solve_cyclic.py

Two things to notice.  First, the **whole** distributed machinery is the single
``solver.solve(communicator=comm)`` line -- everything else is building the system and *checking
the answer*.  Second, the check is real and uses the solver's *own* report rather than a hand-rolled
cutoff: ``report.num_finite_solutions`` is the number of *distinct* finite solutions -- the report
sums :math:`1/\text{multiplicity}` over the finite endpoints, so a genuine **multiple root** (several
paths converging to one true solution of higher multiplicity) is counted once -- and we confirm it
equals the mathematically known cyclic-:math:`n` value.  Just as important, we assert
``report.all_paths_resolved``: had any path failed to track, a genuine root would be silently missing
and the count would come up short, so the report catches the loss rather than letting a quiet 69 pass
for 70.  (A coincident pair of *distinct* paths is the one other way two endpoints meet -- see the
note below.)

.. note::

   Two *distinct* paths landing on the **same** point is a different matter entirely.  Because the
   homotopy's :math:`\gamma` is random, distinct paths meeting is a probability-0 event: if it
   happens it is a **path crossing** -- a sign the tracker under-resolved the paths -- not a benign
   duplicate.  The solver checks for exactly this at the endgame boundary and re-tracks the
   offending paths; with the default predictor and tolerances it essentially never triggers.  The
   :doc:`crossed_paths` tutorial provokes one on purpose and shows what to do about it.

.. _correctness-across-ranks:

.. admonition:: Correctness, not just a path count

   A parallel solve is only useful if it gives the *same answer* as the serial one.  This is
   subtler than it sounds: every rank constructs the homotopy independently, and the start system,
   patch, and :math:`\gamma` are all **random**.  If the ranks disagree on those random choices,
   "path #7" means a different path on each worker and the manager stitches together nonsense --
   the solve finishes, reports a plausible count, and is silently wrong.  Bertini guards against
   this by broadcasting one rank's random seed to all of them at the start of the parallel solve,
   so every rank forms the identical homotopy.  The cyclic check above -- *distinct finite count
   equals the known value* -- is what catches a regression here, which is why one of these two
   examples must verify results and not merely tally paths.

Run the ladder (serial, then 2, 4, 8 workers):

.. code-block:: console

    $ python python/examples/solve_cyclic.py --n 6              # serial baseline
    $ mpirun -n 3 python python/examples/solve_cyclic.py --n 6  # 2 workers
    $ mpirun -n 5 python python/examples/solve_cyclic.py --n 6  # 4 workers
    $ mpirun -n 9 python python/examples/solve_cyclic.py --n 6  # 8 workers

.. BEGIN-TIMING cyclic

.. list-table:: cyclic-6 (720 paths, 156 finite solutions), measured on 8 cores, 2026-06-14 -- run tools/update_scaling_timings.py to refresh
   :header-rows: 1
   :widths: 20 15 25 15

   * - launch
     - workers
     - wall-clock (s)
     - speedup
   * - serial
     - --
     - 357.5
     - 1.0x
   * - ``-n 3``
     - 2
     - 316.6
     - 1.1x
   * - ``-n 5``
     - 4
     - 167.4
     - 2.1x
   * - ``-n 9``
     - 8
     - 120.6
     - 3.0x

.. END-TIMING cyclic

The speedup grows with the worker count but flattens out -- here to about 3x on eight cores, not
the ideal 8x.  That is the *slowest single path* from the opening showing through: once there are
enough workers, the wall-clock is bounded below by the most expensive path, and cyclic-6 at this
tracking tolerance has a few near-singular paths that grind far longer than the rest.  The jump
from 2 to 4 workers is the largest because with only two workers those few heavy paths cannot help
but pile onto one of them; more workers spread them out.  (Note that ``-n 9`` puts nine processes
on eight cores -- the manager only dispatches, so eight workers share the cores with the near-idle
manager.)

Scaling the eigenvalue solve
============================

The eigenvalue solve from :doc:`eigenvalues_by_homotopy` is the opposite regime: a *few* paths,
each *expensive* and of *uneven* cost.  Writing the eigenvector projectively makes the
multihomogeneous Bézout number exactly :math:`n` -- **one path per eigenvalue** -- so an
:math:`n \times n` matrix gives only :math:`n` paths.  To have work worth distributing you scale
the *matrix*, not the path count.

The complete script is :download:`solve_eigenvalues.py <../../../examples/solve_eigenvalues.py>`:

.. literalinclude:: ../../../examples/solve_eigenvalues.py
   :language: python
   :caption: examples/solve_eigenvalues.py

Again the distributed solve is one line, and again it checks correctness -- here against
``numpy.linalg.eigvals``, demanding that every numpy eigenvalue is matched by a homotopy-recovered
one.  Because the matrix is symmetric with integer entries its spectrum is real and (generically)
distinct, so a closest-match comparison is unambiguous.

Run the same ladder on a sizable matrix:

.. code-block:: console

    $ python python/examples/solve_eigenvalues.py --size 24              # serial baseline
    $ mpirun -n 3 python python/examples/solve_eigenvalues.py --size 24  # 2 workers
    $ mpirun -n 5 python python/examples/solve_eigenvalues.py --size 24  # 4 workers
    $ mpirun -n 9 python python/examples/solve_eigenvalues.py --size 24  # 8 workers

.. BEGIN-TIMING eigen

.. list-table:: eigenvalues of a 24x24 symmetric matrix (24 paths), measured on 8 cores, 2026-06-14 -- run tools/update_scaling_timings.py to refresh
   :header-rows: 1
   :widths: 20 15 25 15

   * - launch
     - workers
     - wall-clock (s)
     - speedup
   * - serial
     - --
     - 728.1
     - 1.0x
   * - ``-n 3``
     - 2
     - 477.4
     - 1.5x
   * - ``-n 5``
     - 4
     - 269.0
     - 2.7x
   * - ``-n 9``
     - 8
     - 195.5
     - 3.7x

.. END-TIMING eigen

With only 24 paths the dynamic, demand-driven dispatch earns its keep: under adaptive precision the
per-path cost varies a lot -- a near-collision of eigenvalues makes one path linger at high
precision -- and a *static* split would leave most workers idle while one grinds.  Demand-driven
dispatch keeps the others busy, but it cannot split a *single* path: once enough workers are on
hand, that one slowest eigenvalue sets the floor (here ~195 s, the 3.7x plateau).  Same lesson as
the cyclic case -- distributing helps right up until you hit the longest path, which is exactly the
bound the opening promised.

Threads within ranks: the hybrid model
======================================

Each worker rank can itself track several paths concurrently, one per thread, set with the
``OMP_NUM_THREADS`` environment variable (default 1).  This stacks with MPI: ``R`` worker ranks of
``T`` threads each give ``R × T`` paths in flight.

The launch needs one extra flag, ``--bind-to none``.  By default ``mpirun`` pins each rank to a
single core, which would strangle that rank's threads onto one core; ``--bind-to none`` lets a
rank's threads spread across the machine::

    OMP_NUM_THREADS=4 mpirun -n 3 --bind-to none python python/examples/solve_cyclic.py --n 6
    #               4 threads  x  (3 ranks = 1 manager + 2 workers)  =  8 paths at a time

Why have two knobs for the same cores?  **Reach versus footprint.**  MPI ranks can live on
*different machines* -- that is the only way past one box's core count -- but each rank is a full
process with its own copy of the system.  Threads share one process's memory and stay on one
machine, so they add parallelism without the per-rank memory cost.  The usual recipe: **one rank
per machine (or per NUMA node), threads to fill that machine's cores.**

At a fixed budget of 8 cores, here is the same cyclic-6 solve split between ranks and threads:

.. BEGIN-TIMING hybrid

.. list-table:: cyclic-6 at 8 worker-cores, ranks x threads, measured on 8 cores, 2026-06-14 -- run tools/update_scaling_timings.py to refresh
   :header-rows: 1
   :widths: 40 25 15

   * - layout
     - wall-clock (s)
     - speedup
   * - ``-n 9``, 8 workers x 1 thread
     - 120.6
     - 3.0x
   * - ``-n 5``, 4 workers x 2 threads
     - 126.0
     - 2.8x
   * - ``-n 3``, 2 workers x 4 threads
     - 120.0
     - 3.0x

.. END-TIMING hybrid

On a single shared-memory machine the three layouts land within a few percent of each other -- you
are using the same eight cores either way, and they all bottom out at the same slowest-path floor.
The reason to prefer more *ranks* shows up only across machines (more reach) or under memory
pressure (threads share, ranks duplicate); the reason to prefer more *threads* is to keep the
per-machine process count -- and memory footprint -- low.

Where to go from here
=====================

* The scripts take a ``--n`` / ``--size`` argument: push them up until a single machine is no
  longer enough, then spread the ranks across a cluster with a hostfile
  (``mpirun --hostfile hosts -n ...``).
* Combine with :doc:`parameter_homotopy`: solve once at a generic parameter, then distribute the
  *re-solves* across parameter values -- each is an independent, already-cheap track, ideal for the
  manager-worker pool.
* The same ``solve(communicator=...)`` works for every zero-dim solver class (total degree,
  multihomogeneous, user homotopy), so anything you can solve serially, you can solve distributed.
