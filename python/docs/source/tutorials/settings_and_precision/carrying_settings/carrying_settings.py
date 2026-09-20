"""Carrying settings across related solves -- Bertini 2 tutorial.

Tune one config owner, get_settings(), and set_settings() on each of the rest.
Run:  python carrying_settings.py
"""

import pickle

import numpy as np
import bertini
from bertini import AMPTracker
from bertini.nag_algorithm import ZeroDimConfig
from bertini.tracking import SteppingConfig
from bertini.endgame import EndgameConfig


def build_system():
    """The little system every section solves."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    system = bertini.System()
    system.add_function(x*x + y*y - 1)
    system.add_function(x + y)
    system.add_variable_group(bertini.VariableGroup([x, y]))
    return system


def set_fields_by_name(system):
    """Route fields by name onto whichever config owns them."""
    solver = bertini.ZeroDimSolver(system)
    solver.update(final_tolerance="1e-11",                 # -> EndgameConfig
                  max_num_crossed_path_resolve_attempts=3) # -> ZeroDimConfig

    assert solver.get_config(EndgameConfig).final_tolerance == 1e-11
    assert solver.get_config(ZeroDimConfig).max_num_crossed_path_resolve_attempts == 3

    # Settings the solver keeps in its tracker and its endgame are reached by the same call.
    solver.update(max_step_size="0.05", num_sample_points=6)
    assert solver.get_config(SteppingConfig).max_step_size == bertini.multiprec.real_mp("0.05")
    assert solver.get_config(EndgameConfig).num_sample_points == 6

    # A misspelled field raises immediately rather than doing nothing.
    try:
        solver.update(max_stepsize="0.05")      # no such field anywhere on this solver
        raise AssertionError("should have raised")
    except AttributeError:
        pass

    return solver


def carry_a_whole_bundle(system):
    """get_settings() -> plain picklable dict; set_settings() stamps another solver."""
    reference = bertini.ZeroDimSolver(system)
    reference.update(final_tolerance="1e-11", newton_before_endgame="1e-6")
    settings = reference.get_settings()
    assert set(settings) == set(reference.config_names())     # one entry per config

    # ... later, for each related solve ...
    next_solver = bertini.ZeroDimSolver(system)
    next_solver.set_settings(settings)
    assert next_solver.get_config(EndgameConfig).final_tolerance == 1e-11

    # The bundle is an ordinary picklable value -- store it or ship it to a worker.
    carried = pickle.loads(pickle.dumps(settings))
    worker_solver = bertini.ZeroDimSolver(system)
    worker_solver.set_settings(carried)
    assert worker_solver.get_config(EndgameConfig).final_tolerance == 1e-11

    return settings


def settings_are_precision_agnostic(system, settings):
    """One bundle drops cleanly onto any precision model, and across owner kinds."""
    tuned = bertini.ZeroDimSolver(system, mptype='double')
    tuned.update(final_tolerance="1e-10")
    bundle = tuned.get_settings()

    for mptype in ('multiple', 'adaptive'):
        solver = bertini.ZeroDimSolver(system, mptype=mptype)
        solver.set_settings(bundle)                       # drops on cleanly, any precision model
        assert solver.get_config(EndgameConfig).final_tolerance == 1e-10

    # By default set_settings() applies only the configs the target has and skips the rest.
    tracker = AMPTracker(system)
    tracker.set_settings(settings)             # silently skips the solver-only configs
    try:
        tracker.set_settings(settings, strict=True)
        raise AssertionError("should have raised")
    except KeyError:
        pass


def main():
    system = build_system()
    set_fields_by_name(system)
    settings = carry_a_whole_bundle(system)
    settings_are_precision_agnostic(system, settings)


if __name__ == '__main__':
    main()
