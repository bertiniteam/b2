# This file is part of Bertini 2.
#
# python/test/classes/auxiliary_coordinates_test.py is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
# See the GNU General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""Auxiliary coordinates at the Python layer (b2#403).

A system in null-vector form carries `v` with its own patch, so `v`'s direction is the meaningful
object and its magnitude is an artifact of a normalization nobody chose.  Telling the system that
`v` is auxiliary keeps it out of every judgement -- path truncation, the endgame's divergence
check, finiteness and realness -- without changing anything else about it.

Correctness lives in the C++ tests (``core/test/classes/auxiliary_coordinates_test.cpp`` and the
tracker and solver tests beside them); these check that the surface is reachable and behaves as
documented from Python.
"""

import numpy as np
import pytest

import bertini as pb

_CDT = np.zeros(1, dtype=pb.complex_mp).dtype


def two_groups():
    """(x) and (v): the second standing in for a block whose scale nobody chose."""
    x, v = pb.Variable('x'), pb.Variable('v')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_variable_group(pb.VariableGroup([v]))
    sys.add_function(x * x - 1)
    sys.add_function(v - 1)
    return sys


def point(a, b):
    return np.array([complex(a), complex(b)], dtype=complex)


def test_by_default_every_coordinate_is_judged():
    sys = two_groups()
    assert sys.num_auxiliary_coordinates() == 0
    assert sys.auxiliary_variable_groups() == []
    assert sys.auxiliary_coordinates() == []

    assert sys.is_finite(point(1, 1), 1e5)
    assert not sys.is_finite(point(1, 1e9), 1e5)


def test_an_auxiliary_group_is_not_evidence_about_infinity():
    sys = two_groups()
    sys.set_auxiliary_variable_groups([1])

    assert sys.auxiliary_variable_groups() == [1]
    assert sys.num_auxiliary_coordinates() == 1
    # natural coordinates well inside the working region, the ungoverned block far outside it
    assert sys.is_finite(point(1e3, 1e9), 1e5)
    # and it excludes rather than disables
    assert not sys.is_finite(point(1e9, 1), 1e5)


def test_individual_coordinates_can_be_auxiliary():
    x, v = pb.Variable('x'), pb.Variable('v')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, v]))    # deliberately one group
    sys.add_function(x * x - 1)
    sys.add_function(v - 1)

    sys.set_auxiliary_coordinates([1])
    assert sys.auxiliary_coordinates() == [1]
    assert sys.is_finite(point(1e3, 1e9), 1e5)


def test_realness_is_judged_over_the_same_coordinates():
    sys = two_groups()
    assert not sys.is_real(point(1, complex(1, 0.5)), 1e-8)

    sys.set_auxiliary_variable_groups([1])
    assert sys.is_real(point(1, complex(1, 0.5)), 1e-8)
    assert not sys.is_real(point(complex(1, 0.5), 1), 1e-8)


def test_the_judgements_work_in_multiple_precision():
    sys = two_groups()
    sys.set_auxiliary_variable_groups([1])
    mp = np.array([pb.multiprec.complex_mp('1e3'), pb.multiprec.complex_mp('1e9')], dtype=_CDT)
    assert sys.is_finite(mp, 1e5)


def test_an_index_that_names_nothing_is_refused():
    sys = two_groups()
    with pytest.raises(Exception):
        sys.set_auxiliary_variable_groups([2])
    with pytest.raises(Exception):
        sys.set_auxiliary_coordinates([2])
    assert sys.num_auxiliary_coordinates() == 0


def test_making_every_coordinate_auxiliary_is_refused():
    # a point with nothing to judge would be unconditionally finite and real
    sys = two_groups()
    with pytest.raises(Exception):
        sys.set_auxiliary_variable_groups([0, 1])
    assert sys.num_auxiliary_coordinates() == 0


def test_it_is_part_of_the_system_identity():
    # the tracker truncates on what is NOT auxiliary, so two systems differing only here are
    # different questions and must not recall each other's results
    plain = two_groups()
    with_aux = two_groups()
    with_aux.set_auxiliary_variable_groups([1])
    assert plain.content_digest() != with_aux.content_digest()
