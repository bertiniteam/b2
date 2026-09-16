"""A slice with more linear forms than variables is refused at the boundary (issue #380).

The forms could not be independent, and the orthonormalization used to draw them would hand
back dependent rows without a word.  Interface test: the refusal reaches Python as a
ValueError that names the counts; the C++ tests cover every factory and the boundary cases.
"""

import pytest

import bertini as pb
from bertini import Slice
from bertini.multiprec import complex_mp


def _two_variables():
    x, y = pb.Variable('x'), pb.Variable('y')
    return pb.VariableGroup([x, y])


def test_more_forms_than_variables_is_a_value_error():
    vg = _two_variables()
    with pytest.raises(ValueError, match="more forms than variables"):
        Slice.random_complex(vg, 3)
    with pytest.raises(ValueError, match="more forms than variables"):
        Slice.random_real(vg, 3)
    with pytest.raises(ValueError, match="more forms than variables"):
        Slice.random_complex(vg, 3, True)


def test_a_full_slice_still_works():
    vg = _two_variables()
    assert Slice.random_complex(vg, 2).dimension() == 2
    assert Slice.random_real(vg, 1).dimension() == 1
