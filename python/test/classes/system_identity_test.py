# This file is part of Bertini 2.
#
# python/test/classes/system_identity_test.py is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Bertini 2.
# If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""Interface tests for System content identity (ADR-0042): content_digest / is_same /
seal / is_sealed / intern_system.  Correctness of the identity semantics is tested in C++
(system_identity_test.cpp); these pin the Python surface."""

import pytest
import bertini as pb
from bertini import Variable, VariableGroup


def circle_line():
    x, y = Variable('x'), Variable('y')
    s = pb.System()
    s.add_variable_group(VariableGroup([x, y]))
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x - y)
    return s



def _evaluate_at(system, digits):
    """Materialize `system` at `digits` the only way there is now: evaluate it there.

    There is deliberately no `System.precision(n)` setter -- a system's precision is not a
    property callers manage, it follows the point handed to it (#377).
    """
    import numpy as np
    from bertini.multiprec import complex_mp
    saved = pb.default_precision()
    pb.default_precision(digits)
    try:
        pt = np.array([complex_mp('1', '0')] * system.num_variables())
        system.eval(pt)
    finally:
        pb.default_precision(saved)


def test_content_digest_is_64_hex_and_deterministic():
    a, b = circle_line(), circle_line()
    d = a.content_digest()
    assert len(d) == 64 and int(d, 16) >= 0     # hex, 256 bits
    assert d == b.content_digest()              # independently built, equal content


def test_is_same_and_difference():
    a, b = circle_line(), circle_line()
    assert a.is_same(b) and b.is_same(a)

    x = Variable('x')
    c = pb.System()
    c.add_variable_group(VariableGroup([x]))
    c.add_function(x**3 - 2)
    assert not a.is_same(c)
    assert a.content_digest() != c.content_digest()


def test_seal_blocks_structural_mutation_but_not_eval():
    s = circle_line()
    assert not s.is_sealed()
    s.seal()
    assert s.is_sealed()

    with pytest.raises(Exception):
        s.add_function(Variable('x') + 1)
    with pytest.raises(Exception):
        s.homogenize()

    # transient operations still work
    _evaluate_at(s, 50)
    s.differentiate()
    assert s.is_sealed()


def test_clone_of_sealed_is_unsealed_mutable_copy():
    s = circle_line()
    s.seal()
    c = pb.system.clone(s)
    assert not c.is_sealed()
    assert c.is_same(s)
    c.add_function(Variable('x') - 7)           # no raise: the copy is mutable
    assert not c.is_same(s)


def test_intern_system_unifies_equal_systems():
    rep_a = pb.system.intern_system(circle_line())
    rep_b = pb.system.intern_system(circle_line())
    # one live representative under the hood.  Boost.Python may mint a distinct wrapper
    # object per return, so `is` can't be asserted; shared transient state proves the
    # two handles are one C++ System.
    assert rep_a.is_sealed() and rep_b.is_sealed()
    assert rep_a.is_same(rep_b)
    # Stage a point through ONE handle, then evaluate with no arguments through the OTHER.
    # The no-argument form consumes the system's staged values, so this only works if the two
    # handles share them -- i.e. if they are one C++ System.  (This used to be shown by
    # setting a precision through one handle and reading it back through the other; a System
    # no longer carries a precision, ADR-0057, so the proof now uses the staged point, which
    # is the shared transient state that actually remains.)
    import numpy as np
    rep_a.eval(np.array([complex(3, 0), complex(1, 0)]))
    through_b = [complex(v) for v in np.atleast_1d(np.asarray(rep_b.eval()))]
    assert through_b == [complex(9, 0), complex(2, 0)]

    x = Variable('x')
    other = pb.System()
    other.add_variable_group(VariableGroup([x]))
    other.add_function(x**4 - 3)
    rep_c = pb.system.intern_system(other)
    assert rep_c is not rep_a


def test_digest_ignores_transient_state():
    s = circle_line()
    before = s.content_digest()
    _evaluate_at(s, 60)
    s.differentiate()
    assert s.content_digest() == before


def test_homogenize_changes_digest():
    s = circle_line()
    before = s.content_digest()
    s.homogenize()
    assert s.content_digest() != before
