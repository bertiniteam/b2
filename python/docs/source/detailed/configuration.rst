⚙️ Configuration reference
==========================

Every setting bertini2 has, with its default and what it does. The tables below are read out of
the library when these docs are built, so they describe the version you are reading about rather
than the version someone last wrote a table for.

Three ways to set a setting
---------------------------

Name the field, and it is routed to whichever config holds it. This is the usual way, and it
reaches the solver's own settings, its tracker's, and its endgame's alike::

    solver.update(final_tolerance="1e-11", max_step_size="0.05", num_sample_points=6)

Name the config, when you want to set several of its fields at once::

    solver.configure(stepping={'max_step_size': "0.05"},
                     endgame={'final_tolerance': "1e-8"})

Hand over a whole config object, which is what :meth:`get_settings` and :meth:`set_settings` carry
between solvers::

    settings = tuned_solver.get_settings()
    next_solver.set_settings(settings)

A misspelled name raises immediately and suggests the one you meant, rather than doing nothing.

Numbers, exactly
----------------

Numeric settings take strings, and a string is the noise-free way to write one: ``"1e-11"`` means
exactly that, where the Python float ``1e-11`` is the nearest double to it. For the settings
stored as exact rationals -- a sample ladder's ratio, the homotopy times -- a Python float is
refused outright, since ``0.1`` the double is not one tenth; write ``"1/10"`` or ``"0.1"`` and
both give the exact rational, or pass an ``int``, a :class:`fractions.Fraction`, or a
:class:`bertini.multiprec.rational_mp`.

.. _settings-without-a-config:

Tracker settings that have no config struct
-------------------------------------------

These are set by calling a method on the tracker rather than by naming a field, so they appear in
no config listing, including the generated ones below. They are among the most consequential
settings there are, and this is the only place they are collected.

.. list-table::
   :header-rows: 1
   :widths: 38 20 42

   * - Setting
     - Default
     - What it does
   * - ``tracker.infinite_truncation(bool)``
     - ``True``
     - Whether to truncate such paths at all.
   * - ``tracker.precision_preservation(bool)``
     - ``False``
     - Whether an adaptive tracker returns to its starting precision when a path ends.
   * - ``tracker.reinitialize_initial_step_size(bool)``
     - ``True``
     - Whether each path starts from the configured initial step size rather than inheriting the
       last one used.
   * - ``tracker.set_max_wall_clock_duration(seconds)``
     - none
     - A wall-clock deadline; a track in progress is abandoned between steps once it passes.
       Cleared with ``clear_max_wall_clock_time()``. See
       :doc:`/tutorials/settings_and_precision/stopping_and_budgets/index`.

The configs
-----------

Everything below is generated. Each heading is a config class; the keyword beside it is what
:meth:`configure` takes, and every field name is what :meth:`update` accepts directly.

.. config-reference::
