🌟 Everyday API
****************

The handful of names you reach for constantly, grouped by what you are doing.  Everything here lives
at the **top level** of the package -- ``import bertini`` and tab-complete ``bertini.`` -- unless a
submodule is named.  For the exhaustive, auto-generated reference (every class, every method, every
submodule), see :doc:`the full API index <api>`.

.. note::

   Casing is a signal.  A **lowercase** ``_mp`` name is a concrete arbitrary-precision *number*
   (like a numpy dtype); a **CapWords** name in ``bertini.symbolics`` is a *symbol* you build
   expressions from.  So ``bertini.complex_mp('0.1')`` is a value, while ``bertini.symbolics.Complex``
   is a coefficient node.


Building a system
=================

Start from variables, combine them with arithmetic and the operators below, and collect functions
into a :class:`System`.

``Variable('x')``
    A single symbolic variable.
``variables('x', 5)``
    A list of indexed variables ``x0..x4`` (pass a range for other index sets).
``VariableGroup([x, y])`` / ``VariableGroup('x', 5)``
    An (affine) group of variables -- pass a list, or a name and a count.
``System()``
    The polynomial system.  Add to it with ``sys.add_function(f)``, ``sys.add(grp, f, g)``, and the
    block builders ``sys.add_functions(...)``, ``sys.add_linear(...)``,
    ``sys.add_products_of_linears(...)``, ``sys.randomize(...)``.
``coefficient(v)`` / ``coefficients(A)``
    Turn exact values (ints, ``fractions.Fraction``, exact strings, multiprec numbers) into
    coefficient nodes -- refuses Python floats, which would cap precision.
``jacobian([f, g], [x, y])``
    The symbolic Jacobian, as a numpy object array of expressions.
``random_matrix(m, n)``
    A random (numeric or symbolic) coefficient matrix -- e.g. a random linear projection.

See the tutorials :doc:`/tutorials/zerodim_solver/index` and
:doc:`/tutorials/user_product_of_linears/index`.


Symbols, constants, and operators
==================================

The math vocabulary for *building* expressions on variables.  These are at the top level; the
:mod:`bertini.operators` module re-exports just this vocabulary so you can
``from bertini.operators import *`` without pulling in the rest of the package.

``E``, ``Pi``, ``I``
    The symbolic constants (``I`` is the imaginary unit).
``sin cos tan asin acos atan exp log sqrt``
    Elementary functions that build expression nodes, e.g. ``sin(x) + Pi*y``.
``bertini.symbolics``
    The full, flat symbolic namespace -- ``symbolics.Variable``, ``symbolics.Complex``,
    ``symbolics.NamedExpression``, the operator node types (``symbolics.Sum``, ``symbolics.Power``,
    ...), and ``symbolics.AbstractNode`` for isinstance-based tree walking.  You rarely construct
    these by hand -- literals auto-convert -- but reach here when you want to *know* you are holding
    a symbol.
``Named(expr, 'a')``
    Give a subexpression a name (a ``NamedExpression``).


Numbers and precision
=====================

Arbitrary-precision numbers, usable directly and as numpy dtypes.  Also in
:mod:`bertini.multiprec` (which additionally carries the numeric ``sin``/``cos``/... that act on
numbers rather than symbols).

``complex_mp`` / ``real_mp`` / ``int_mp`` / ``rational_mp``
    Arbitrary-precision complex / real / integer / rational number types.
``default_precision(n)``
    Get or set the global working precision (decimal digits).  **Global mutable state** -- set it
    before building the values whose precision you care about.

See :doc:`/tutorials/precision_matters/index` and :doc:`/tutorials/precision_models/index`.


Solving
=======

The zero-dimensional solve and the homotopy building blocks.

``ZeroDimSolver(sys, ...)``
    Solve a square system for its isolated solutions; ``.solve()`` then ``.solutions()`` /
    ``.all_solutions()``.  The big everyday entry point.
``HomotopySolver`` / ``SolutionPathCollector``
    Track a user-supplied homotopy, and collect solution paths.
``Slice`` / ``Slice.from_coefficients(coeffs, vars)``
    A linear slice (the linear part of a witness set); build one from an exact coefficient matrix.
``StartSystemType``
    Which start system to use (e.g. total degree).  The homotopy helpers
    (``parameter_sweep``, ``moving_homotopy``, ``coefficient_parameter_homotopy``,
    ``blend_homotopy``) stay in :mod:`bertini.nag_algorithm`.

See :doc:`/tutorials/zerodim_solver/index`, :doc:`/tutorials/all_solutions/index`, and
:doc:`/tutorials/parameter_homotopy/index`.


Tracking and endgames
======================

Lower-level path tracking, when you want to drive it yourself.

``AMPTracker`` / ``DoublePrecisionTracker`` / ``MultiplePrecisionTracker``
    The path trackers.  ``AMPTracker`` (adaptive multiple precision) is the usual choice.
``Predictor``
    Predictor method enum (``Euler``, ``HeunEuler``, ``RK4``, ...).
``SuccessCode``
    The result enum every track / solve step returns (``Success``, ...).
``bertini.endgame``
    The endgames for singular endpoints: ``AMPCauchyEndgame``, ``AMPPowerSeriesEndgame``, and the
    fixed-precision variants (``FixedDouble...``, ``FixedMultiple...``).

See :doc:`/tutorials/tracking_nonsingular/index` and :doc:`/tutorials/manual_endgame_usage/index`.


Settings
========

Every tracker/solver owns its configuration structs (in :mod:`bertini.tracking`,
:mod:`bertini.endgame`, :mod:`bertini.nag_algorithm`).  Set individual fields by name in one call:

``owner.set(**fields)`` / ``owner.update(**fields)``
    Route each named setting to whichever config owns it, e.g.
    ``solver.set(final_tolerance='1e-11')`` or
    ``tracker.get_stepping().set(max_step_size='0.05')``.  ``set`` and ``update`` are the same.

See :doc:`/tutorials/carrying_settings/index`.


Other namespaces
================

``bertini.parse``
    Read classic Bertini-1 input files / strings into a :class:`System`.
``bertini.parallel``
    MPI helpers for distributed solves.
``bertini.random``
    Seeding (``set_random_seed``) and the random draws behind ``random_matrix``.
``bertini.logging``
    Logging configuration.

Enums live at the root for convenience: ``SuccessCode``, ``Predictor``, ``MonomialOrder``,
``StartSystemType`` (they also remain in their submodules).
