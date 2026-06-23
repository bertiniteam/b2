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
from bertini.tracking import SteppingConfig, NewtonConfig
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


def test_update_accepts_string_for_mpfr_field():
    # strings are the blessed noise-free input for numeric settings
    c = SteppingConfig().update(max_step_size="0.05")
    assert c.max_step_size == pb.multiprec.Float("0.05")


def test_update_accepts_mpfr_float_for_mpfr_field():
    c = SteppingConfig().update(max_step_size=pb.multiprec.Float("0.25"))
    assert c.max_step_size == pb.multiprec.Float("0.25")


def test_update_rejects_python_float_for_mpfr_field():
    # floats stay rejected by policy: 0.05 the double is not 1/20
    with pytest.raises(TypeError):
        SteppingConfig().update(max_step_size=0.05)


def test_update_accepts_string_for_double_tolerance_field():
    # NumErrorT (tolerance) fields are plain doubles -- a high-precision number would have no value
    # there -- but update() still takes the same noise-free string spelling as the exact fields, so
    # the whole config surface is uniform.  (Regression: these used to throw Boost ArgumentError.)
    assert TolerancesConfig().update(final_tolerance="1e-11").final_tolerance == 1e-11
    assert TolerancesConfig().update(newton_before_endgame="1e-7").newton_before_endgame == 1e-7


def test_update_accepts_string_for_min_step_size():
    # min_step_size is a plain double (def_readwrite) where its siblings are mpq_rational; update()
    # papers over that difference so a string works on every stepping field alike.
    from bertini.tracking import SteppingConfig as SC
    assert SC().update(min_step_size="1e-50").min_step_size == 1e-50


def test_update_accepts_string_for_integer_field():
    from bertini.tracking import NewtonConfig
    assert NewtonConfig().update(max_num_newton_iterations="4").max_num_newton_iterations == 4


def test_update_float_still_works_on_double_fields():
    # the float-rejection policy is only for the exact (mpq/mpfr) fields; a plain double is the
    # natural input for a NumErrorT tolerance and must keep working.
    assert TolerancesConfig().update(final_tolerance=1e-11).final_tolerance == 1e-11


def test_zero_dim_config_is_precision_agnostic():
    # ZeroDimConfig is no longer templated on the complex type: there is one config, keyed 'zero_dim'
    # on every precision model (was zero_dim_config_double_prec / _multiprec).  This is what lets a
    # carried settings bundle apply across a double, multiple, or adaptive solver unchanged.
    from bertini.nag_algorithm import ZeroDimConfig
    import bertini.nag_algorithm as na
    assert not hasattr(na, 'ZeroDimConfigDoublePrec')
    assert not hasattr(na, 'ZeroDimConfigMultiprec')

    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_function(x * x + y * y - 1); s.add_function(x + y)
    s.add_variable_group(pb.VariableGroup([x, y]))
    for mptype in ('double', 'multiple', 'adaptive'):
        names = pb.nag_algorithm.ZeroDim(s, mptype=mptype).config_names()
        assert 'zero_dim' in names, names


def test_zero_dim_config_times_accept_strings():
    # the homotopy times are stored precision-free (mpq_rational) but exposed as real and take the
    # same string spelling as every other numeric field.
    from bertini.nag_algorithm import ZeroDimConfig
    c = ZeroDimConfig().update(endgame_boundary="0.05", start_time="1", target_time="0")
    assert c.endgame_boundary == pb.multiprec.Float("0.05")


def test_regeneration_slice_tolerances_are_prefixed():
    # RegenerationConfig's tolerances are the slice-MOVING tracking tolerances (Bertini 1's SliceTol*
    # family), distinct from the main tracking tolerances in TolerancesConfig.  They carry a slice_
    # prefix so every config field name is unique across structs -- the precondition for routing a
    # field to its config without naming the struct.
    from bertini.nag_algorithm import RegenerationConfig
    c = RegenerationConfig().update(slice_newton_before_endgame="1e-7",
                                    slice_newton_during_endgame="1e-8",
                                    slice_final_tolerance="1e-12")
    assert c.slice_newton_before_endgame == 1e-7
    assert c.slice_final_tolerance == 1e-12
    # the un-prefixed names belong only to TolerancesConfig now
    assert not hasattr(c, 'newton_before_endgame')
    assert not hasattr(c, 'final_tolerance')


def test_no_field_name_collisions_across_configs():
    # The slice_ rename leaves every config field name unique across all of an owner's configs, which
    # is what lets a field be routed to its owning config unambiguously.
    from bertini.nag_algorithm import (ZeroDimCauchyAdaptivePrecisionTotalDegree as ZD,
                                       TolerancesConfig, RegenerationConfig)
    from bertini.config import writable_fields
    tol = set(writable_fields(TolerancesConfig))
    regen = set(writable_fields(RegenerationConfig))
    assert tol & regen == set(), "tolerances/regeneration still share field names: {}".format(tol & regen)


def test_update_rejects_garbage_string():
    with pytest.raises(Exception):
        SteppingConfig().update(max_step_size="not a number")


def test_string_value_roundtrips_through_to_dict():
    c = SteppingConfig().update(max_step_size="0.125")
    d = c.to_dict()
    assert d['max_step_size'] == pb.multiprec.Float("0.125")
    assert SteppingConfig.from_dict(d) == c


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


# ------------------------------------------------ owner.update(**fields) field router

def test_owner_update_routes_fields_by_name(solver):
    # the headline ergonomic: set fields on the owner without naming the config struct -- each field
    # goes to whichever config owns it.
    from bertini.nag_algorithm import ZeroDimConfig
    solver.update(final_tolerance="1e-11", max_num_crossed_path_resolve_attempts=3)
    assert solver.get_config(TolerancesConfig).final_tolerance == 1e-11
    assert solver.get_config(ZeroDimConfig).max_num_crossed_path_resolve_attempts == 3


def test_owner_update_is_chainable(solver):
    assert solver.update(final_tolerance="1e-9") is solver


def test_owner_update_accepts_strings(solver):
    # strings work through the router for every numeric field, same as the per-config update().
    solver.update(final_tolerance="1e-12")
    assert solver.get_config(TolerancesConfig).final_tolerance == 1e-12


def test_tracker_update_routes_to_its_configs(tracker):
    from bertini.tracking import SteppingConfig, NewtonConfig
    tracker.update(max_step_size="0.05", max_num_newton_iterations=2)
    assert tracker.get_config(SteppingConfig).max_step_size == pb.multiprec.Float("0.05")
    assert tracker.get_config(NewtonConfig).max_num_newton_iterations == 2


def test_owner_update_routes_only_across_this_owners_configs(solver):
    # max_step_size lives on the tracker's SteppingConfig, not on the algorithm's configs -- the
    # router refuses it on the algorithm rather than silently doing nothing.
    with pytest.raises(AttributeError):
        solver.update(max_step_size="0.05")


def test_owner_update_rejects_unknown_field(solver):
    with pytest.raises(AttributeError):
        solver.update(finaltol="1e-9")


# ------------------------------------------------ get_settings / set_settings (the carry)

def _square():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System()
    s.add_function(x ** 2 + y ** 2 - 1); s.add_function(x + y)
    s.add_variable_group(pb.VariableGroup([x, y]))
    return s


def test_get_settings_is_a_named_dict_of_configs():
    a = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    settings = a.get_settings()
    assert set(settings) == set(a.config_names())
    assert isinstance(settings['tolerances'], TolerancesConfig)


def test_settings_round_trip_onto_another_solver():
    from bertini.nag_algorithm import ZeroDimConfig
    a = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    a.update(final_tolerance="1e-12", max_num_crossed_path_resolve_attempts=4)

    b = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    b.set_settings(a.get_settings())
    assert b.get_config(TolerancesConfig).final_tolerance == 1e-12
    assert b.get_config(ZeroDimConfig).max_num_crossed_path_resolve_attempts == 4


def test_settings_carry_across_precision_models():
    # the de-templated, precision-agnostic configs are what make this work: a bundle from a multiple-
    # precision solver applies unchanged to a double or adaptive one.  This is the cross-stage carry
    # an NID-style workflow needs.
    src = pb.nag_algorithm.ZeroDim(_square(), mptype='multiple')
    src.update(final_tolerance="1e-11")
    for mptype in ('double', 'adaptive'):
        dst = pb.nag_algorithm.ZeroDim(_square(), mptype=mptype)
        dst.set_settings(src.get_settings())
        assert dst.get_config(TolerancesConfig).final_tolerance == 1e-11


def test_settings_bundle_is_picklable():
    import pickle
    a = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    a.update(final_tolerance="1e-9")
    restored = pickle.loads(pickle.dumps(a.get_settings()))
    b = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    b.set_settings(restored)
    assert b.get_config(TolerancesConfig).final_tolerance == 1e-9


def test_set_settings_skips_inapplicable_by_default_strict_raises():
    from bertini.tracking import AMPTracker
    settings = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square()).get_settings()
    trk = AMPTracker(_square())            # a tracker has no 'tolerances' / 'zero_dim'
    trk.set_settings(settings)             # non-strict: silently skips them
    with pytest.raises(KeyError):
        trk.set_settings(settings, strict=True)


def test_set_settings_accepts_dict_of_fields():
    a = ZeroDimCauchyAdaptivePrecisionTotalDegree(_square())
    a.set_settings({'tolerances': {'final_tolerance': '1e-10'}})
    assert a.get_config(TolerancesConfig).final_tolerance == 1e-10
