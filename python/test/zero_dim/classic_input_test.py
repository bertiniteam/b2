"""solver.to_classic_input(system): emit a Bertini 1 classic input file whose CONFIG
reflects THIS solver's tracking settings, and whose INPUT is the natural system you pass.

The system-level System.to_classic_input(**kwargs) (knobs supplied by hand, for sweeps) is
tested in classes/classic_writer_test.py; round-trip fidelity is pinned in C++
(classic_parsing_test::classic_writer_round_trips_a_system).  Here we check that the *algorithm*
supplies the CONFIG -- in particular the predictor and precision mode, which a bare System
cannot know.
"""

import bertini as pb
from bertini import Predictor


def _circle_line():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return sys


def test_emits_config_and_input_sections():
    sys = _circle_line()
    text = pb.ZeroDimSolver(sys).to_classic_input(sys)
    assert 'CONFIG' in text
    assert 'INPUT' in text
    assert text.count('END;') == 2          # one closes CONFIG, one closes INPUT
    assert 'tracktype: 0;' in text          # zero-dim solve


def test_passed_system_supplies_input():
    sys = _circle_line()
    text = pb.ZeroDimSolver(sys).to_classic_input(sys)
    assert 'variable_group x, y;' in text
    assert 'function f0, f1;' in text
    assert 'f0 = ' in text and 'f1 = ' in text


def test_default_predictor_is_rkf45():
    # the precision-mode default is exercised in test_precision_mode_follows_the_solver
    # (and differs by branch); the predictor default is RKF45 regardless.
    sys = _circle_line()
    text = pb.ZeroDimSolver(sys).to_classic_input(sys)
    assert 'odepredictor: 5;' in text       # RKF45, the Bertini 2 default


def test_precision_mode_follows_the_solver():
    sys = _circle_line()
    assert 'mptype: 0;' in pb.ZeroDimSolver(sys, mptype='double').to_classic_input(sys)
    assert 'mptype: 1;' in pb.ZeroDimSolver(sys, mptype='multiple').to_classic_input(sys)
    assert 'mptype: 2;' in pb.ZeroDimSolver(sys, mptype='adaptive').to_classic_input(sys)


def test_predictor_change_is_reflected():
    """The discriminator: a bare System cannot know the predictor; only the solver's tracker can.
    Change it and the emitted odepredictor must follow."""
    sys = _circle_line()
    solver = pb.ZeroDimSolver(sys)
    solver.get_tracker().predictor(Predictor.Euler)
    text = solver.to_classic_input(sys)
    assert 'odepredictor: 0;' in text       # Euler
    assert 'odepredictor: 5;' not in text


def test_tolerances_come_through():
    sys = _circle_line()
    text = pb.ZeroDimSolver(sys).to_classic_input(sys)
    # Bertini 2 defaults: newton-before 1e-5, newton-during 1e-6, final 1e-11.
    assert 'tracktolbeforeeg: 1e-05;' in text
    assert 'tracktolduringeg: 1e-06;' in text
    assert 'finaltol: 1e-11;' in text
