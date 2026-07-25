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


def test_float_coefficients_print_full_precision():
    """A coefficient must round-trip EXACTLY through classic input: streaming at the
    ostream default (6 significant digits) silently truncated every printed system
    (found via a cellcap bundle whose crit-curve system was a 1e-6 impostor)."""
    import bertini
    bertini.default_precision(30)
    third = bertini.multiprec.real_mp(1) / bertini.multiprec.real_mp(3)
    x, = bertini.variables(['x'])
    s = bertini.System()
    s.add_variable_group([x])
    s.add_function(x - bertini.coefficient(third))
    txt = s.to_classic_input()
    assert '0.333333333333333333333333333' in txt, txt   # full digits, not 0.333333
    import numpy as np
    s2 = bertini.parse.system(txt)
    # the reparsed coefficient is the IDENTICAL binary value
    val = s2.eval(np.array([bertini.multiprec.complex_mp(third)]))
    assert float(abs(complex(np.asarray(val).ravel()[0]))) == 0.0

    # complex coefficients: both components at full precision through the pair form
    c = bertini.multiprec.complex_mp(third, -third)
    s3 = bertini.System()
    s3.add_variable_group([x])
    s3.add_function(x - bertini.coefficient(c))
    s4 = bertini.parse.system(s3.to_classic_input())
    val = s4.eval(np.array([c]))
    assert float(abs(complex(np.asarray(val).ravel()[0]))) == 0.0
