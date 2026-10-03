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

clone is a memory-isolating copy (C++ ``Clone``, ADR-0027): it SHARES the immutable node DAG
(variables included -- variables are canonical by name anyway) and only gives the copy its own
evaluation memory, so a clone is safe to evaluate concurrently with the original yet still shares
its variables.  Sharing variables is what makes clone-then-concatenate work (issue #256): the two
systems have, by construction, identical variable orderings.  These tests check that a clone (and a
concatenation) is not just structurally the right shape but *evaluates* identically -- functions AND
Jacobian, with and without a path variable -- and that adding functions to one system does not leak
into the other.  Systems are also pickleable (a boost-archive round-trip), so copy.copy /
copy.deepcopy / pickle round-trip a System too; those are exercised here alongside clone.
"""

import copy

import numpy as np
import pytest

import bertini as pb


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
    s.add_linear(np.array([[2, 1], [1, -3]]), np.array([x, y]), [-1, 4])

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


def test_concatenate_of_cloned_system_issue_256():
    # Regression for issue #256: clone a system to reuse its variable-group setup, give the
    # original and the clone different functions, then concatenate.  The old serialization-based
    # clone minted fresh variable nodes, so the two systems' orderings compared unequal and
    # concatenate threw "differing variable orderings".  The ADR-0027 clone shares the variable
    # nodes, so the orderings match and concatenate succeeds (and evaluates correctly).
    x, y = pb.Variable('x'), pb.Variable('y')
    orig = pb.System()
    orig.add_variable_group(pb.VariableGroup([x, y]))

    new = pb.system.clone(orig)              # reuse the variable-group setup, no re-declaration
    orig.add_function(2 * x - 7 * y - 1)
    new.add_function(x - 3 * y)

    u = pb.system.concatenate(new, orig)     # <- issue #256 threw here
    assert u.num_functions() == 2

    pt = np.array([complex(2.0), complex(5.0)])   # x=2, y=5
    # row 0 = new's  x - 3y    = -13 ; row 1 = orig's 2x - 7y - 1 = -32
    _assert_close(_eval(u, pt), [complex(-13.0), complex(-32.0)])


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


@pytest.mark.parametrize("roundtrip", [
    lambda s: pickle.loads(pickle.dumps(s)),
    copy.copy,
    copy.deepcopy,
    pb.system.clone,
])
def test_a_copy_combines_with_its_original(roundtrip):
    # Loading from an archive re-interns (ADR-0068), so a pickled, copied or deep-copied System is
    # over THE SAME variables as its original and concatenates with it.  It used to come back over
    # fresh variables of the same names, and was refused: "differing variable orderings".
    a, b = _two_systems()
    u = pb.system.concatenate(a, roundtrip(b))
    assert u.num_functions() == 3

    pt = np.array([complex(2.0, 0.3), complex(-1.0, 0.5)])
    _assert_close(_eval(u, pt), _eval(a, pt) + _eval(b, pt))
    assert [str(v) for v in u.variable_ordering()] == ['x', 'y']


def test_separately_homogenized_deep_copies_build_a_homotopy():
    # The cellular port's projective move: three deep copies, homogenized one by one, the fixed one
    # patched.  straight_line_homotopy assembles .target and .start by concatenation, and refused.
    x, y = pb.Variable('x'), pb.Variable('y')

    def over_xy(f):
        s = pb.System()
        s.add_variable_group(pb.VariableGroup([x, y]))
        s.add_function(f)
        return s

    fixed = copy.deepcopy(over_xy(x * x + y * y - 1))
    start = copy.deepcopy(over_xy(x - 2))
    end = copy.deepcopy(over_xy(y - 3))
    fixed.homogenize()
    fixed.auto_patch()
    start.homogenize()
    end.homogenize()

    built = pb.nag_algorithm.straight_line_homotopy(end, start, fixed=fixed, gamma=1)
    assert built.homotopy.num_functions() == built.target.num_functions()
    assert [str(v) for v in built.target.variable_ordering()] == \
           [str(v) for v in fixed.variable_ordering()]


def _over(variables, function):
    s = pb.System()
    s.add_variable_group(pb.VariableGroup(list(variables)))
    s.add_function(function)
    return s


def test_a_different_order_is_refused_everywhere():
    # The homotopy builders used to count variables, and blended (x, y) with (y, x) by position:
    # y - 3 at (x, y) = (7, 11) evaluated to 4.
    x, y = pb.Variable('x'), pb.Variable('y')
    xy = _over([x, y], x - 2)
    yx = _over([y, x], y - 3)
    fixed = _over([x, y], x + y - 5)

    for combine in (lambda: pb.system.concatenate(xy, yx),
                    lambda: pb.nag_algorithm.straight_line_homotopy(yx, xy, projectivize=False),
                    lambda: pb.nag_algorithm.straight_line_homotopy(yx, xy, fixed=fixed,
                                                                    projectivize=False)):
        with pytest.raises(RuntimeError, match='same variables in a different order'):
            combine()


def test_the_same_order_blends_the_rows_it_was_given():
    x, y = pb.Variable('x'), pb.Variable('y')
    fixed = _over([x, y], x + y - 5)
    built = pb.nag_algorithm.straight_line_homotopy(_over([x, y], y - 3), _over([x, y], x - 2),
                                                    fixed=fixed, gamma=1, projectivize=False)
    at_the_end = _eval(built.homotopy, np.array([complex(7.0), complex(11.0)]), complex(0.0))
    _assert_close(at_the_end, [complex(13.0), complex(8.0)])       # x + y - 5, and y - 3


def test_a_different_grouping_and_different_variables_are_refused():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    one_group = _over([x, y, z], x - 1)
    two_groups = pb.System()
    two_groups.add_variable_group(pb.VariableGroup([x, y]))
    two_groups.add_variable_group(pb.VariableGroup([z]))
    two_groups.add_function(z - 1)
    with pytest.raises(RuntimeError, match='the grouping does not'):
        pb.system.concatenate(one_group, two_groups)
    with pytest.raises(RuntimeError, match='different variables'):
        pb.system.concatenate(_over([x, y], x - 1), _over([x, z], z - 1))


def test_a_system_that_declares_no_variables_takes_the_others():
    x, y = pb.Variable('x'), pb.Variable('y')
    declared = _over([x, y], x * x + y * y - 1)
    bare = pb.System()
    bare.add_function(x - y)
    assert bare.num_variables() == 0

    pt = np.array([complex(2.0), complex(5.0)])
    _assert_close(_eval(pb.system.concatenate(declared, bare), pt), [complex(28.0), complex(-3.0)])
    _assert_close(_eval(pb.system.concatenate(bare, declared), pt), [complex(-3.0), complex(28.0)])


# ----------------------------- nodes compare by identity -----------------------------

def test_a_variable_made_twice_is_one_variable():
    # Variables are canonical by name, so both handles are on one node; Python's default ==
    # compared the two wrapper objects and said they differed.
    assert pb.Variable('x') == pb.Variable('x')
    assert not (pb.Variable('x') != pb.Variable('x'))
    assert pb.Variable('x') != pb.Variable('y')
    assert hash(pb.Variable('x')) == hash(pb.Variable('x'))


def test_equal_expressions_are_one_node():
    x, y = pb.Variable('x'), pb.Variable('y')
    assert x * y + 1 == 1 + y * x
    assert x * y != x + y
    assert len({x, pb.Variable('x'), y}) == 2
    assert {x: 'first'}[pb.Variable('x')] == 'first'


def test_a_node_is_not_equal_to_something_that_is_not_a_node():
    x = pb.Variable('x')
    assert not (x == 1)
    assert x != 1
    assert not (x == 'x')
    assert not (x == None)      # noqa: E711 -- the comparison is the thing under test


def test_a_variable_group_compares_by_its_variables():
    x, y = pb.Variable('x'), pb.Variable('y')
    assert pb.VariableGroup([x, y]) == [pb.Variable('x'), pb.Variable('y')]
    assert pb.VariableGroup([x, y]) != [y, x]
    assert _over([x, y], x - 1).variable_ordering() == [x, y]


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
