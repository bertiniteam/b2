"""The classic-input round trip, and the three defects found in it.

``System.to_classic_input()`` emits a complete Bertini 1 file (``CONFIG ... END;`` then
``INPUT ... END;``) because its purpose is running the same problem in Bertini 1.
``bertini.parse.system()`` reads the INPUT-section body.  The two must round-trip, and when
they cannot, they must say so rather than hand back an empty system.

Covers bertiniteam/b2 #395, #396 and #397.
"""

import re

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.symbolics import Variable


def whitney_ish():
    """f = (x^2+y^2+z^2+3)^2 - 16(x^2+y^2)  -- a torus.  TOTAL DEGREE 4."""
    x, y, z = Variable('x'), Variable('y'), Variable('z')
    s = System()
    vg = VariableGroup()
    for v in (x, y, z):
        vg.append(v)
    s.add_variable_group(vg)
    s.add_functions([(x**2 + y**2 + z**2 + 3)**2 - 16 * (x**2 + y**2)])
    return s, vg


# --------------------------------------------------------------- #395

def test_empty_system_emits_classic_input_without_crashing():
    """REGRESSION (#395): this SEGFAULTED.

    ``CoefficientBound`` ended in ``max(f_vals...maxCoeff(), dh_dx...maxCoeff(), bound)``,
    and with no functions both are EMPTY -- Eigen's ``maxCoeff()`` on an empty array is
    undefined behaviour.  ``to_classic_input()`` reaches it because the CONFIG section
    emits ``coefficientbound:``.
    """
    text = System().to_classic_input()
    assert 'coefficientbound:' in text


def test_system_with_variables_but_no_functions_emits_classic_input():
    """The same crash, with variables present: the trigger is zero FUNCTIONS."""
    s = System()
    vg = VariableGroup()
    vg.append(Variable('x'))
    s.add_variable_group(vg)
    assert 'coefficientbound:' in s.to_classic_input()


# --------------------------------------------------------------- #396

def test_full_classic_file_round_trips():
    """REGRESSION (#396): parse.system() could not read what to_classic_input() writes.

    It failed at line 2 on the CONFIG block: ``expected "=" found "tracktype: 0;"``.
    """
    s, _ = whitney_ish()
    back = bertini.parse.system(s.to_classic_input())
    assert back.num_functions() == 1
    assert back.num_variables() == 3


def test_round_trip_preserves_the_input_section_verbatim():
    """The CONFIG section is NOT compared: ``coefficientbound`` is estimated by evaluating
    at random points, so it differs between two emissions of the very same system."""
    s, _ = whitney_ish()
    a = s.to_classic_input()
    b = bertini.parse.system(a).to_classic_input()
    body = lambda t: t.split('INPUT', 1)[1]
    assert body(a) == body(b)


def test_bare_input_section_body_still_parses():
    """The pre-existing calling convention must keep working -- the unwrapping is a no-op
    for text that is already a bare body."""
    p = bertini.parse.system("variable_group x, y, z;\nfunction f0;\nf0 = x^2+y^2+z^2-1;\n")
    assert p.num_functions() == 1
    assert p.num_variables() == 3


def test_input_wrapper_alone_is_accepted():
    p = bertini.parse.system(
        "INPUT\nvariable_group x, y, z;\nfunction f0;\nf0 = x^2+y^2+z^2-1;\nEND;\n")
    assert p.num_functions() == 1
    assert p.num_variables() == 3


@pytest.mark.parametrize("bad", [
    "this is not a system at all",
    "",
    "INPUT\nEND;\n",
])
def test_unparseable_input_raises_rather_than_returning_an_empty_system(bad):
    """REGRESSION (#396), the dangerous half.

    The binding discarded ``parse()``'s return value, and ``parse()`` leaves its result
    untouched when it matches nothing -- so Python received a structurally valid, entirely
    EMPTY System and no error.  A caller round-tripping through text carried on with zero
    functions.  Silence is the bug; an exception is the fix.
    """
    with pytest.raises(RuntimeError):
        bertini.parse.system(bad)


# --------------------------------------------------------------- #397

@pytest.mark.parametrize("expr,total_degree", [
    ("x^4", 4),
    ("(x+y)^2", 2),
    ("(x^2+y^2)^2", 4),
    ("(x^2+y^2+z^2+3)^2", 4),
])
def test_parsed_power_of_a_sum_reports_its_total_degree(expr, total_degree):
    """REGRESSION (#397): PowerOperator::Degree(VariableGroup) SUMMED the per-variable
    degrees, which is valid only for a monomial.  ``(x+y)^2`` has degree 2 in x and 2 in y
    but TOTAL degree 2; it reported 4, and ``(x^2+y^2+z^2+3)^2`` reported 12 instead of 4.

    Invisible from Python's ``**`` (which builds an IntegerPowerOperator, always correct)
    and invisible in ungrouped ``degrees()`` -- it showed only through
    ``Degrees(Variables())``, which is what ``DegreeBound()`` calls, and DegreeBound feeds
    AMP.
    """
    p = bertini.parse.system(
        "variable_group x, y, z;\nfunction f0;\nf0 = %s;\n" % expr)
    assert int(re.search(r'degreebound: (\d+)', p.to_classic_input()).group(1)) == total_degree


def test_degree_bound_agrees_between_a_built_and_a_parsed_system():
    """The consequence that matters: a system round-tripped through classic input must
    track under the SAME AMP degree bound as the one it was built from."""
    s, _ = whitney_ish()
    p = bertini.parse.system(s.to_classic_input())
    db = lambda S: int(re.search(r'degreebound: (\d+)', S.to_classic_input()).group(1))
    assert db(s) == db(p) == 4
