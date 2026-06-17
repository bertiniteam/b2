# This file is part of Bertini 2.
#
# python/test/classes/clone_concatenate_test.py is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team

"""Tests for bertini.system.clone and bertini.system.concatenate.

clone is a *serialization round-trip* deep copy (boost text-archive, see C++ ``Clone``), so it is
the natural place for serialization bugs to surface -- the same round-trip the distributed solver
uses to broadcast systems.  These tests check that a clone (and a concatenation) is not just
structurally the right shape but *evaluates* identically -- functions AND Jacobian, with and without
a path variable -- and is a genuinely independent deep copy.  Systems are now pickleable (a
boost-archive round-trip, same machinery as clone), so copy.copy / copy.deepcopy / pickle round-trip
a System too; those are exercised here alongside clone.
"""

import copy

import numpy as np
import pytest

import bertini as pb
from bertini import linalg


TOL = 1e-12


def _eval(s, pt, time=None):
    return [complex(c) for c in (s.eval(pt) if time is None else s.eval(pt, time))]


def _jac(s, pt, time=None):
    s.differentiate()
    j = s.eval_jacobian(pt) if time is None else s.eval_jacobian(pt, time)
    # eval_jacobian returns a 1-D gradient for a single-function system and a 2-D matrix otherwise;
    # ravel() flattens both.  Clone and original share a shape, so a flat compare is valid.
    return [complex(z) for z in np.asarray(j).ravel()]


def _assert_close(a, b, tol=TOL):
    assert len(a) == len(b)
    for u, v in zip(a, b):
        assert abs(u - v) < tol, (u, v)


# ----------------------------- clone: evaluation parity -----------------------------

def test_clone_preserves_function_values():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x + y * y - 1)
    s.add_function(x * y - 2)

    c = pb.system.clone(s)
    assert c.num_functions() == 2 and c.num_variables() == 2

    pt = np.array([complex(2.3, -1.1), complex(-0.7, 0.4)])
    _assert_close(_eval(c, pt), _eval(s, pt))


def test_clone_preserves_jacobian():
    # The Jacobian is the fragile part of the round trip: C++ Clone re-derives derivatives
    # (Differentiate()) because the serialized SLP's derivative outputs did not survive faithfully.
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x * y)
    s.add_function(x - y * y)

    c = pb.system.clone(s)
    pt = np.array([complex(1.5, 0.3), complex(-2.0, 1.0)])
    _assert_close(_jac(c, pt), _jac(s, pt))


def test_clone_with_path_variable_preserves_eval_and_jacobian():
    # A time-dependent system exercises the time-derivative path the C++ Clone comment specifically
    # flags as having read stale memory after the round trip.
    x, y, t = pb.Variable('x'), pb.Variable('y'), pb.Variable('t')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_path_variable(t)
    s.add_function((1 - t) * (x * x + y * y - 1) + t * (x - y))

    c = pb.system.clone(s)
    assert c.have_path_variable()

    sp = np.array([complex(2.0, 0.5), complex(3.0, -0.25)])
    tm = complex(0.4, -0.2)
    _assert_close(_eval(c, sp, tm), _eval(s, sp, tm))
    _assert_close(_jac(c, sp, tm), _jac(s, sp, tm))


def test_clone_of_linear_forms_block_preserves_eval_and_jacobian():
    # A LinearFormsBlock (add_linear) is a structured evaluation block; it must survive the
    # serialization round trip just like the products-of-linears block does.
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x + y * y - 1)
    linalg.add_linear(s, np.array([[2, 1], [1, -3]]), np.array([x, y]), [-1, 4])

    c = pb.system.clone(s)
    assert list(c.degrees()) == list(s.degrees())

    pt = np.array([complex(0.6, 0.2), complex(-0.9, 0.1)])
    _assert_close(_eval(c, pt), _eval(s, pt))
    _assert_close(_jac(c, pt), _jac(s, pt))


def test_clone_of_homogenized_patched_preserves_eval():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x + y * y - 1)
    s.homogenize()
    s.auto_patch()

    c = pb.system.clone(s)
    assert c.num_functions() == s.num_functions()
    assert c.num_variables() == s.num_variables()

    # one coordinate per (now 3) variables, on the patch is not required for a plain eval check
    pt = np.array([complex(0.3, 0.1), complex(-0.5, 0.2), complex(0.8, -0.4)])
    _assert_close(_eval(c, pt), _eval(s, pt))


# ----------------------------- clone: independence (deep copy) -----------------------------

def test_clone_is_independent_of_the_original():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x + y * y - 1)

    c = pb.system.clone(s)
    pt = np.array([complex(2.0), complex(3.0)])
    before = _eval(c, pt)

    # mutate the original after cloning
    s.add_function(x - y)
    assert s.num_functions() == 2
    assert c.num_functions() == 1, "clone must not see the original's new function"
    _assert_close(_eval(c, pt), before)


# ----------------------------- concatenate -----------------------------

def _two_systems():
    x, y = pb.Variable('x'), pb.Variable('y')
    a = pb.System()
    a.add_variable_group(pb.VariableGroup([x, y]))
    a.add_function(x * x + y * y - 1)
    b = pb.System()
    b.add_variable_group(pb.VariableGroup([x, y]))
    b.add_function(x - y)
    b.add_function(x * y - 2)
    return a, b


def test_concatenate_appends_and_evaluates():
    a, b = _two_systems()
    u = pb.system.concatenate(a, b)
    assert u.num_functions() == a.num_functions() + b.num_functions() == 3

    pt = np.array([complex(2.0, 0.3), complex(-1.0, 0.5)])
    _assert_close(_eval(u, pt), _eval(a, pt) + _eval(b, pt))


def test_concatenate_is_independent_of_inputs():
    a, b = _two_systems()
    u = pb.system.concatenate(a, b)
    pt = np.array([complex(1.0), complex(2.0)])
    before = _eval(u, pt)

    x = pb.Variable('x')
    a.add_function(x)        # mutate an input after concatenating
    assert u.num_functions() == 3, "concatenation must not see the input's new function"
    _assert_close(_eval(u, pt), before)


# ----------------------------- pickle / copy / deepcopy round-trip -----------------------------

import pickle


def _time_system():
    x, y, t = pb.Variable('x'), pb.Variable('y'), pb.Variable('t')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_path_variable(t)
    s.add_function((1 - t) * (x * x + y * y - 1) + t * (x - y))
    return s


@pytest.mark.parametrize("roundtrip", [
    lambda s: pickle.loads(pickle.dumps(s)),
    copy.copy,
    copy.deepcopy,
])
def test_pickle_copy_deepcopy_preserve_eval_and_jacobian(roundtrip):
    # System is pickleable (boost-archive round-trip), so copy.copy / copy.deepcopy / pickle all
    # produce an independent System that evaluates identically -- functions and Jacobian, including
    # the time-derivative path that the C++ Clone comment flags as fragile.
    s = _time_system()
    c = roundtrip(s)

    sp = np.array([complex(2.0, 0.5), complex(3.0, -0.25)])
    tm = complex(0.4, -0.2)
    _assert_close(_eval(c, sp, tm), _eval(s, sp, tm))
    _assert_close(_jac(c, sp, tm), _jac(s, sp, tm))


def test_pickle_roundtrip_is_independent_of_the_original():
    # A pickled-and-restored System must be a genuine deep copy, not an alias.
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_variable_group(pb.VariableGroup([x, y]))
    s.add_function(x * x + y * y - 1)

    c = pickle.loads(pickle.dumps(s))
    pt = np.array([complex(2.0), complex(3.0)])
    before = _eval(c, pt)

    s.add_function(x - y)
    assert c.num_functions() == 1, "restored System must not see the original's new function"
    _assert_close(_eval(c, pt), before)
