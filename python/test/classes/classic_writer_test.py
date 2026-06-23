"""System.to_classic_input(): emit a Bertini 1 classic input file (CONFIG + INPUT).

Round-trip fidelity (emit -> re-parse -> same evaluations) is pinned in C++
(classic_parsing_test::classic_writer_round_trips_a_system).  Here we check the binding: the
wrapper, the system body, and the selectable CONFIG knobs.
"""

import bertini


def _circle_line():
    x, y = bertini.Variable('x'), bertini.Variable('y')
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x * x + y * y - 1)
    sys.add_function(x - y)
    return sys


def test_emits_config_and_input_sections():
    text = _circle_line().to_classic_input()
    assert 'CONFIG' in text
    assert 'INPUT' in text
    assert text.count('END;') == 2                 # one closes CONFIG, one closes INPUT
    assert 'tracktype: 0;' in text                 # zero-dim solve


def test_emits_the_system_body():
    text = _circle_line().to_classic_input()
    assert 'variable_group x, y;' in text
    assert 'function f0, f1;' in text
    assert 'f0 = ' in text and 'f1 = ' in text


def test_mptype_and_predictor_are_selectable():
    assert 'mptype: 2;' in _circle_line().to_classic_input()             # adaptive default
    assert 'mptype: 0;' in _circle_line().to_classic_input(mptype=0)     # double
    assert 'odepredictor: 0;' in _circle_line().to_classic_input(odepredictor=0)
    assert 'odepredictor: 5;' in _circle_line().to_classic_input()       # RKF45 default
