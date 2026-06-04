# This file is part of Bertini 2.
#
# python/test/nid/nid_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/nid/nid_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/nid/nid_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#  silviana amethyst
#  UWEC
#

"""Tests for the NumericalIrreducibleDecomposition framework/skeleton: that the
algorithm class is exported, constructs, is wired to its configs via the reusable
ConfiguredVisitor interface (no NID-specific Python code), exposes a tracker/endgame
and a result datatype, and that its (not-yet-implemented) compute entry point raises."""

import unittest

import bertini as pb
from bertini.nag_algorithm import (
    NIDCauchyAdaptivePrecision,
    NIDPowerSeriesDoublePrecision,
    RegenerationConfig,
    TolerancesConfig,
    NumericalIrreducibleDecompositionMultiplePrecision,
    WitnessSetMultiplePrecision,
)


def _square_system():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x + y)
    s.add_variable_group(pb.VariableGroup([x, y]))
    return s


class NIDFrameworkTest(unittest.TestCase):
    def setUp(self):
        self.nid = NIDCauchyAdaptivePrecision(_square_system())

    def test_constructs_and_is_exported(self):
        self.assertIsNotNone(self.nid)
        # a second concrete instantiation should exist too
        self.assertIsNotNone(NIDPowerSeriesDoublePrecision(_square_system()))

    def test_exposes_its_configs_via_the_reusable_interface(self):
        names = self.nid.config_names()
        self.assertIn('regeneration', names)
        self.assertIn('tolerances', names)
        self.assertIn('sharpening', names)
        self.assertIn('post_processing', names)

    def test_set_and_get_config_roundtrips(self):
        self.nid.set_config(RegenerationConfig().update(start_level=3))
        self.assertEqual(self.nid.get_config(RegenerationConfig).start_level, 3)

    def test_configure_many_subconfigs_in_one_call(self):
        self.nid.configure(regeneration={'start_level': 2},
                           tolerances={'final_tolerance': 1e-11})
        self.assertEqual(self.nid.get_config(RegenerationConfig).start_level, 2)
        self.assertEqual(self.nid.get_config(TolerancesConfig).final_tolerance, 1e-11)

    def test_has_tracker_and_endgame(self):
        self.assertIsNotNone(self.nid.get_tracker())
        self.assertIsNotNone(self.nid.get_endgame())

    def test_solve_not_yet_implemented(self):
        with self.assertRaises(RuntimeError):
            self.nid.solve()

    def test_default_decomposition_is_empty(self):
        decomp = self.nid.decomposition()
        self.assertEqual(decomp.num_witness_sets(), 0)
        self.assertEqual(list(decomp.nonempty_codimensions()), [])


class NIDDataTypeTest(unittest.TestCase):
    def test_result_and_witness_set_construct(self):
        decomp = NumericalIrreducibleDecompositionMultiplePrecision()
        self.assertEqual(decomp.num_witness_sets(), 0)

        ws = WitnessSetMultiplePrecision()
        self.assertEqual(ws.degree(), 0)


if __name__ == '__main__':
    unittest.main()
