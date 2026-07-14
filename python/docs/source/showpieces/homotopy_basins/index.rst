🌀 Homotopy Basins
******************

**Every pixel of this image is a zero-dimensional solve.**  Rasterize a window of the complex
coefficient plane of :math:`f(x;c) = x^9 - 9x - c`; at each pixel, run the genuine total-degree
gamma-trick homotopy

.. math::

   H(x,t) \;=\; \gamma\, t\,(x^9 - 1) \;+\; (1-t)\,(x^9 - 9x - c)

and track all nine start roots from :math:`t=1` to :math:`t=0`.  Nothing below is a synthetic
texture — every channel is tracked data.

.. image:: homotopy_basins.png
   :width: 100%
   :align: center
   :alt: neon comet cluster: eight glowing discriminant arcs with rainbow spiral vortices at their
         curled roots, striped basins streaming to the right

The **hue** is the phase of a landing *fingerprint* :math:`\sum_k \zeta^k x_k` over the nine
destinations: smooth inside a basin, jumping exactly where the homotopy's swept discriminant
permutes which root each path reaches.  The **brightness** is the tracker's own step count — the
blazing arcs are the set of :math:`c` for which the straight-line homotopy passes through a
singular system at some :math:`t`, and the adaptive stepper piles up tiny steps exactly there.
The **stripes** are level bands of the fingerprint's magnitude, and the rainbow whirlpools are
its zeros, coiled at the roots of the arcs.

The gamma trick, photographed
=============================

Why arcs, and why do they all stream the same way?  The family's branch points — the
:math:`c`-values where two roots of :math:`f` collide — form the ring
:math:`c = -8\zeta,\ \zeta^8 = 1`.  The homotopy drags every pixel's system along a straight line
through *system space*, and that line grazes the discriminant along one arc per branch point:

.. math::

   c(\mu) \;=\; -\,\frac{4\,\zeta}{(1+\mu)^{1/4}} \;-\; \mu,
   \qquad \mu \in \gamma\cdot(0,\infty),

each arc leaving its branch point and trailing to infinity in the :math:`-\gamma` direction.
That is the **gamma trick** made visible: the arcs are precisely the *bad* target systems for
this choice of :math:`\gamma`, a measure-zero set that a different :math:`\gamma` re-aims
elsewhere — rotate :math:`\gamma` and the whole comet cluster swings around the ring.  For a
generic pixel the paths pass *near* but never *through* the discriminant, and the closer the
brush, the harder the tracker works: the glow you see is `num_total_steps_taken` doing numerical
algebraic geometry the hard way.  The rare white speckle marks pixels where a path genuinely
failed — the gamma trick's measure-zero fine print, caught on camera.

Start simple: five roots, four arcs
===================================

The same construction for :math:`f(x;c) = x^5 - 5x - c` keeps every feature legible: four branch
points at :math:`c = -4\zeta,\ \zeta^4=1`, one glowing arc each, five basins:

.. image:: homotopy_basins_teaching.png
   :width: 85%
   :align: center
   :alt: the degree-5 version: four glowing arcs with spiral vortices, five striped basins

Crossing an arc means the straight-line homotopy passed the discriminant on one side rather than
the other, so the destinations of exactly two of the tracked paths trade places — the hue field
jumps by the corresponding transposition's fingerprint.  Walking a small loop *around* an arc's
endpoint (the branch point itself) is precisely the monodromy loop of
:doc:`the Monodromy Loom </showpieces/monodromy_loom/index>` — the two showpieces are the same
mathematics seen from parameter space and from solution space.

How it is built
===============

One exact homotopy per pixel (pixel coordinates are snapped to rationals — every coefficient the
function tree sees is exact, and there is no randomness anywhere, so the tracked data is fully
deterministic):

.. literalinclude:: homotopy_basins.py
   :language: python
   :start-at: def _track_pixel(
   :end-at: return (lands

The rest — the fingerprint reduction, the stress normalization, the banding and bloom — is
rendering; see ``homotopy_basins.py`` in full.

.. note::

   This is a raster showpiece: PNG only (the image *is* a raster of solves; there is no
   meaningful vector form), regenerated through ``tools/refresh_doc_artifacts.py``.  It is the
   heaviest doc artifact — about 17 million tracked paths across the two frames (roughly 10–15
   minutes on 12 cores).  It is not a doctest — run the generator directly with
   ``python homotopy_basins.py``.
