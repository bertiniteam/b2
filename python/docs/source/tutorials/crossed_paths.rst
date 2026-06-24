🪢 Crossed paths: when two paths become one
*********************************************

Homotopy continuation tracks one path per start point, from a generic start system at :math:`t=1`
to your target system at :math:`t=0`.  The homotopy carries a random parameter :math:`\gamma`, and
that randomness is load-bearing: it guarantees, *with probability 1*, that the paths stay
**separate** for all :math:`0 < t \le 1`.  Distinct paths simply do not meet.

So if two distinct paths *do* arrive at the **same** point, something has gone wrong.  It is not a
benign "duplicate" -- it is a **path crossing**: somewhere along the way the numerical tracker took
a step coarse enough to jump from one path onto a neighbouring one, and from there the two are
indistinguishable.  They reach the endgame boundary as the same point, the endgame sends both to
the same solution, and you are left with a duplicate where there should have been two answers --
**a solution silently lost.**

This tutorial does three things: provoke a crossing on purpose, see how Bertini 2 detects and
*fixes* it, and -- the real point -- show why you should never see one with the defaults.

What Bertini does at the endgame boundary
=========================================

A zero-dim solve has a definite shape: track every path to the **endgame boundary** (by default
:math:`t = 0.1`), then run the endgame from there to :math:`t=0`.  In between sits the **midpath
check**.  It compares every pair of boundary points and, if two agree to within a relaxed
same-point tolerance, flags a crossing.  Then :func:`EGBoundaryAction` re-tracks the offending
paths with **tightened settings** -- a smaller tolerance *and* a higher-order predictor -- and
checks again, up to ``max_num_crossed_path_resolve_attempts`` times (default 2).

Two distinct paths agreeing to several digits at :math:`t=0.1` is a measure-zero coincidence, which
is exactly why it is a reliable alarm: under correct tracking it never fires.

Provoking one
=============

The classic way to under-resolve paths is a low-order predictor with a loose tolerance.  We use the
**Euler** predictor (order 1) and a loose pre-endgame tolerance on the cyclic-5 system (120 paths,
70 distinct finite solutions).  The complete script is :download:`crossed_paths.py
<../../../examples/crossed_paths.py>`:

.. literalinclude:: ../../../examples/crossed_paths.py
   :language: python
   :caption: examples/crossed_paths.py

The key knobs are three lines on the solver:

.. code-block:: python

    solver.get_tracker().predictor(pb.tracking.Predictor.Euler)   # order-1 predictor: the culprit

    tols = solver.get_config(pb.nag_algorithm.TolerancesConfig)
    tols.newton_before_endgame = 1e-4                              # loose tracking before the endgame
    solver.set_config(tols)

    zd = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    zd.max_num_crossed_path_resolve_attempts = 0                  # 0 = detect & report, do NOT re-track
    solver.set_config(zd)

.. note::

   The fixed seed only makes the *bad behaviour* reproducible so we can talk about it.  Pinning a
   seed is **never** the fix for a flaky solve -- the fixes are tolerance, precision, and predictor.
   See :doc:`solving_at_scale`.

Seeing the damage, then the repair
===================================

Run with re-tracking **disabled** first, so the crossing stands uncorrected:

.. code-block:: text

    Euler + loose tolerance, re-tracking DISABLED (max_num_crossed_path_resolve_attempts=0):
      crossings detected : 2
      crossed path indices: [16, 58]
      midpath check passed: False
      distinct finite solutions: 69 / 70

      --> the crossed paths merged into one, so a true solution was lost

Paths 16 and 58 jumped onto each other.  The solver *tells you* this -- the report from
:func:`endgame_boundary_metadata` carries ``passed = False`` and the offending indices -- but with
re-tracking off it does nothing about it, and the distinct finite count comes up one short: **69
instead of 70.**

Now the default behaviour, re-tracking **enabled**:

.. code-block:: text

    Same solve, re-tracking ENABLED (the default):
      crossings detected : 2
      re-track attempts  : 1
      midpath check passed: True
      distinct finite solutions: 70 / 70

      --> the crossing was resolved and the lost solution recovered.

One re-track attempt -- tighter tolerance, and Euler bumped up to the default RKF45 -- pulls the two
paths back apart, and the missing solution reappears.  **70 of 70.**

Reading the report
==================

Every solve exposes what the midpath check found, via :func:`endgame_boundary_metadata`:

.. code-block:: python

    report = solver.endgame_boundary_metadata()
    report.passed                  # did the check ultimately pass (no crossings left)?
    report.num_crossings_detected  # how many crossed paths were found on the first check
    report.num_resolve_attempts    # how many re-track attempts were performed
    report.crossed_path_indices    # which paths were involved

Use it as a correctness gate in your own scripts: if ``report.passed`` is ``False`` after a solve,
some crossings were left unresolved and the affected solutions may be wrong -- re-solve with a
better predictor or a tighter tolerance.

What to actually do about crossings
===================================

* **Prefer prevention.**  The defaults exist precisely so this never happens: the default predictor
  is **RKF45** (order 4 with an embedded error estimate), not Euler.  Swap the one Euler line out --
  or just don't put it in -- and the crossing disappears at the source.  A modestly tighter
  ``newton_before_endgame`` does the same.
* **Let the solver self-correct.**  With ``max_num_crossed_path_resolve_attempts > 0`` (the default
  is 2) the re-track usually fixes a stray crossing for you, as above.
* **Decide what "give up" means.**  Set ``max_num_crossed_path_resolve_attempts = 0`` to *detect and
  report only* -- the solve keeps tracking the crossed paths through the endgame (you never get
  fewer answers than you would without the check), but it does not spend time re-tracking; you read
  ``report.passed`` and decide.  This is the conservative choice when you would rather re-solve the
  whole system with better settings than retry individual paths.

The honest summary: a crossing is a signal you tracked too coarsely.  Bertini 2 will catch it and,
within a bounded number of attempts, usually fix it -- but the best solve is the one where the
predictor and tolerance are good enough that the alarm never sounds.
