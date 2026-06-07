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

import pytest

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


# ---------------------------------------------- pure config-struct conveniences

def test_update_is_chainable_and_sets_fields():
    c = SteppingConfig().update(min_num_steps=3, max_num_steps=999)
    assert isinstance(c, SteppingConfig)
    assert c.min_num_steps == 3
    assert c.max_num_steps == 999


def test_update_rejects_unknown_field():
    with pytest.raises(AttributeError):
        SteppingConfig().update(not_a_real_field=5)


def test_repr_shows_class_and_fields():
    text = repr(SteppingConfig().update(min_num_steps=7))
    assert text.startswith("SteppingConfig(")
    assert "min_num_steps=7" in text


def test_to_dict_and_from_dict_roundtrip():
    c = NewtonConfig().update(max_num_newton_iterations=4, min_num_newton_iterations=2)
    d = c.to_dict()
    assert d['max_num_newton_iterations'] == 4
    assert NewtonConfig.from_dict(d) == c


def test_equality_by_value():
    a = SteppingConfig().update(min_num_steps=5)
    b = SteppingConfig().update(min_num_steps=5)
    c = SteppingConfig().update(min_num_steps=6)
    assert a == b
    assert a != c


# -------------------------------------------------------------- tracker configs

@pytest.fixture
def tracker():
    return AMPTracker(_square_system())


def test_config_names_lists_the_typelist(tracker):
    names = tracker.config_names()
    assert 'stepping' in names
    assert 'newton' in names
    # the AMP tracker's precision config is AdaptiveMultiplePrecisionConfig (AMPConfig)
    assert 'amp' in names


def test_set_config_then_get_config_roundtrips(tracker):
    tracker.set_config(NewtonConfig().update(max_num_newton_iterations=9))
    assert tracker.get_config(NewtonConfig).max_num_newton_iterations == 9


def test_configure_one_call_many_subconfigs(tracker):
    tracker.configure(stepping={'min_num_steps': 11},
                      newton={'max_num_newton_iterations': 6})
    assert tracker.get_config(SteppingConfig).min_num_steps == 11
    assert tracker.get_config(NewtonConfig).max_num_newton_iterations == 6


def test_get_stepping_internal_ref_updates_in_place(tracker):
    # the legacy mutable accessor + update() edits the live config
    tracker.get_stepping().update(min_num_steps=21)
    assert tracker.get_config(SteppingConfig).min_num_steps == 21


# ------------------------------------------------------------ algorithm configs
# Previously there was no way to touch an algorithm's configs from Python.

@pytest.fixture
def solver():
    return ZeroDimCauchyAdaptivePrecisionTotalDegree(_square_system())


def test_algorithm_exposes_its_configs(solver):
    names = solver.config_names()
    assert 'tolerances' in names
    assert 'post_processing' in names


def test_set_and_get_algorithm_config(solver):
    tol = solver.get_config(TolerancesConfig).update(final_tolerance=1e-11)
    solver.set_config(tol)
    assert solver.get_config(TolerancesConfig) == tol
