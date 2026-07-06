"""Intro to slices -- witness sets, slices, and start systems in Bertini 2.

Run:  python intro_to_slices.py
"""

import pickle

import numpy as np

import bertini
from bertini import Slice
from bertini import multiprec as mp
from bertini.nag_algorithm import WitnessSetMultiplePrecision


def build_slices():
    """The slice: a stack of linear forms."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    # two linear forms on (x, y):  2x + 3y + 1   and   x - y + 4
    s = Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])

    assert s.dimension() == 2        # two forms -> cuts a 2-dimensional component
    assert s.num_variables() == 2

    # or generate a generic one (the usual case -- a witness set's slice is random)
    x0, x1, x2 = bertini.Variable('x0'), bertini.Variable('x1'), bertini.Variable('x2')
    vg = bertini.VariableGroup([x0, x1, x2])

    generic = Slice.random_complex(vg, 2)    # two random hyperplanes in 3-space
    assert generic.dimension() == 2
    assert generic.num_variables() == 3

    return x, y, s, x0, x1, x2, vg, generic


def sequence_semantics(x, y, s):
    """A slice is a sequence of forms -- and the shape rules matter."""
    # an ELEMENT: slice[i] is the i-th form's coefficient VECTOR (1-D), "a line is a vector"
    v = np.asarray(s[0])
    assert v.shape == (3,)
    assert [complex(c).real for c in v] == [2.0, 3.0, 1.0]      # 2x + 3y + 1

    # iterating yields the form vectors
    assert len(list(s)) == 2

    # a SUB-COLLECTION: slice[i:j] is a new Slice
    first = s[:1]
    assert first.dimension() == 1

    # the whole matrix: always 2-D, (num_forms, num_variables + 1)
    assert s.coefficients().shape == (2, 3)

    # a single-form slice never collapses to 1-D
    one_form = Slice.from_coefficients([[2, 3, 1]], [x, y])
    assert one_form.coefficients().shape == (1, 3)     # NOT (3,) -- it never collapses
    assert np.asarray(one_form[0]).shape == (3,)       # the vector view is explicit


def slice_rides_along(s, x0, x1, x2, vg, generic):
    """A slice rides along with a system."""
    # as_system(): a standalone System of just the slice's forms
    only_slice = s.as_system()
    assert only_slice.num_functions() == 2

    # add_to(): append the slice's forms to an existing system (over the same variables)
    sphere = bertini.System()
    sphere.add_variable_group(vg)
    sphere.add_function(x0*x0 + x1*x1 + x2*x2 - 1)      # the unit sphere: a surface in 3-space
    generic.add_to(sphere)
    assert sphere.num_functions() == 3                  # 1 sphere equation + 2 slice forms


def build_witness_set(x0, x1, x2, vg):
    """The witness set: assemble the triple."""
    def pt(*entries):
        return np.array([mp.complex_mp(str(e)) for e in entries], dtype=mp.complex_mp)

    sys = bertini.System()
    sys.add_variable_group(vg)
    sys.add_function(x0*x0 + x1*x1 + x2*x2 - 1)         # the sphere again (a 2-dim component)
    slice2 = Slice.random_complex(vg, 2)

    # all at once: points + slice + system
    w = WitnessSetMultiplePrecision([pt(1, 0, 0), pt(0, 1, 0)], slice2, sys)
    assert w.degree() == 2          # two witness points
    assert w.dimension() == 2       # a surface
    assert w.is_consistent()        # 3 variables - 1 equation == slice dimension 2

    # or incrementally
    w2 = WitnessSetMultiplePrecision()
    w2.set_system(sys)
    w2.set_slice(slice2)
    w2.add_point(pt(1, 0, 0))
    w2.add_point(pt(0, 1, 0))
    assert w2.degree() == 2

    # a witness set prints a readable summary
    print(repr(w))

    return w


def serialize_and_emit(w):
    """Carrying it around, and emitting to Bertini 1."""
    w_again = pickle.loads(pickle.dumps(w))
    assert w_again.degree() == 2
    assert w_again.is_consistent()

    # emit the witness (square) system as a Bertini 1 classic input file
    square = w.witness_system()          # system + slice, the square system to track
    assert square.num_functions() == 3
    text = w.to_classic_input()
    assert 'CONFIG' in text and 'INPUT' in text


def start_from_slices(x, y):
    """Building start systems from slices (regeneration)."""
    start = bertini.System()
    start.add_variable_group(bertini.VariableGroup([x, y]))
    sa = Slice.from_coefficients([[1, 0, -1], [1, 0, 1]], [x, y])    # (x - 1)(x + 1)
    sb = Slice.from_coefficients([[0, 1, -1], [0, 1, -2]], [x, y])   # (y - 1)(y - 2)
    start.add_slices_as_products([sa, sb])
    assert list(start.degrees()) == [2, 2]


def main():
    x, y, s, x0, x1, x2, vg, generic = build_slices()
    sequence_semantics(x, y, s)
    slice_rides_along(s, x0, x1, x2, vg, generic)
    w = build_witness_set(x0, x1, x2, vg)
    serialize_and_emit(w)
    start_from_slices(x, y)


if __name__ == '__main__':
    main()
