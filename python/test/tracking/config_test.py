# This file is part of Bertini 2.
#
# python/test/tracking/config_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/tracking/config_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/tracking/config_test.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""Tests for the Pythonic config ergonomics (bertini.config) on trackers and
nag_algorithms: update()/repr/to_dict/from_dict/eq on configs, and
config_types()/get_config()/set_config()/configure() on owners."""

import unittest

import bertini as pb
from bertini.tracking import AMPTracker
from bertini.tracking.config import SteppingConfig, NewtonConfig
from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionTotalDegree, TolerancesConfig


def _square_system():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_function(x**2 + y**2 - 1)
    s.add_function(x + y)
    s.add_variable_group(pb.VariableGroup([x, y]))
    return s


class ConfigStructTest(unittest.TestCase):
    """The pure config-struct conveniences, independent of any owner."""

    def test_update_is_chainable_and_sets_fields(self):
        c = SteppingConfig().update(min_num_steps=3, max_num_steps=999)
        self.assertIsInstance(c, SteppingConfig)
        self.assertEqual(c.min_num_steps, 3)
        self.assertEqual(c.max_num_steps, 999)

    def test_update_rejects_unknown_field(self):
        with self.assertRaises(AttributeError):
            SteppingConfig().update(not_a_real_field=5)

    def test_repr_shows_class_and_fields(self):
        text = repr(SteppingConfig().update(min_num_steps=7))
        self.assertTrue(text.startswith("SteppingConfig("))
        self.assertIn("min_num_steps=7", text)

    def test_to_dict_and_from_dict_roundtrip(self):
        c = NewtonConfig().update(max_num_newton_iterations=4, min_num_newton_iterations=2)
        d = c.to_dict()
        self.assertEqual(d['max_num_newton_iterations'], 4)
        self.assertEqual(NewtonConfig.from_dict(d), c)

    def test_equality_by_value(self):
        a = SteppingConfig().update(min_num_steps=5)
        b = SteppingConfig().update(min_num_steps=5)
        c = SteppingConfig().update(min_num_steps=6)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)


class TrackerConfigTest(unittest.TestCase):
    def setUp(self):
        self.tracker = AMPTracker(_square_system())

    def test_config_names_lists_the_typelist(self):
        names = self.tracker.config_names()
        self.assertIn('stepping', names)
        self.assertIn('newton', names)
        # the AMP tracker's precision config is AdaptiveMultiplePrecisionConfig (AMPConfig)
        self.assertIn('amp', names)

    def test_set_config_then_get_config_roundtrips(self):
        self.tracker.set_config(NewtonConfig().update(max_num_newton_iterations=9))
        self.assertEqual(self.tracker.get_config(NewtonConfig).max_num_newton_iterations, 9)

    def test_configure_one_call_many_subconfigs(self):
        self.tracker.configure(stepping={'min_num_steps': 11},
                               newton={'max_num_newton_iterations': 6})
        self.assertEqual(self.tracker.get_config(SteppingConfig).min_num_steps, 11)
        self.assertEqual(self.tracker.get_config(NewtonConfig).max_num_newton_iterations, 6)

    def test_get_stepping_internal_ref_updates_in_place(self):
        # the legacy mutable accessor + update() edits the live config
        self.tracker.get_stepping().update(min_num_steps=21)
        self.assertEqual(self.tracker.get_config(SteppingConfig).min_num_steps, 21)


class AlgorithmConfigTest(unittest.TestCase):
    """Previously there was no way to touch an algorithm's configs from Python."""

    def setUp(self):
        self.solver = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square_system())

    def test_algorithm_exposes_its_configs(self):
        names = self.solver.config_names()
        self.assertIn('tolerances', names)
        self.assertIn('post_processing', names)

    def test_set_and_get_algorithm_config(self):
        tol = self.solver.get_config(TolerancesConfig).update(final_tolerance=1e-11)
        self.solver.set_config(tol)
        self.assertEqual(self.solver.get_config(TolerancesConfig), tol)


if __name__ == '__main__':
    unittest.main()
