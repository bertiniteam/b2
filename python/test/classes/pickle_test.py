# This file is part of Bertini 2.
#
# python/test/classes/pickle_test.py is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team

"""Pickling support across the bound types.

What is pickleable, and how:

  * ``multiprec.real_mp`` / ``multiprec.complex_mp`` -- C++ boost-archive pickle suite (exact, every bit
    of the mantissa survives -- see the precision-sensitive round-trip below);
  * config structs (stepping, newton, tolerances, the zero-dim configs, ...) -- via their
    ``to_dict``/``from_dict`` enhancement (``config.py``);
  * Boost.Python enums (``SuccessCode``, ``Predictor``, ...) -- via a ``copyreg`` registration that
    reconstructs the member from its int value (their built-in ``__reduce_ex__`` is broken);
  * result metadata (``SolutionMetaData*``) -- enhanced like a config, and now round-trips because
    its enum fields (e.g. ``endgame_success_code``) are pickleable.

(System pickling lives in ``clone_concatenate_test.py``, next to the clone serialization tests.)
"""

import copy
import pickle

import pytest

import bertini as pb
import bertini.multiprec as mp


# ----------------------------- numbers -----------------------------

@pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)
def test_complex_pickles_exactly_at_high_precision(precision):
    # The whole point of multiprecision: a float pickle that drops low bits is silently wrong.  Use
    # a value with more significant digits than double can hold and demand exact equality.
    z = mp.complex_mp("1.23456789012345678901234567890123456789",
                   "-9.87654321098765432109876543210987654321")
    assert pickle.loads(pickle.dumps(z)) == z
    assert copy.deepcopy(z) == z


@pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)
def test_float_pickles_exactly_at_high_precision(precision):
    f = mp.real_mp("3.14159265358979323846264338327950288419716939937510")
    assert pickle.loads(pickle.dumps(f)) == f
    assert copy.deepcopy(f) == f


# ----------------------------- enums -----------------------------

@pytest.mark.parametrize("member", [
    pb.SuccessCode.Success,
    pb.SuccessCode.GoingToInfinity,
    pb.Predictor.RK4,
    pb.Predictor.Euler,
])
def test_enum_members_pickle(member):
    r = pickle.loads(pickle.dumps(member))
    assert r == member
    assert int(r) == int(member)
    assert type(r) is type(member)


def test_enum_is_not_config_enhanced():
    # Boost.Python enums subclass int; they must NOT be mistaken for config structs (that corrupts
    # their repr and pickling).  The clean enum repr is the tell.
    SC = pb.SuccessCode
    assert not getattr(SC, "_b2_config_enhanced", False)
    assert "SuccessCode.Success" in repr(SC.Success)


# ----------------------------- config structs -----------------------------

@pytest.mark.parametrize("cls", [
    pb.tracking.SteppingConfig,
    pb.tracking.NewtonConfig,
    pb.tracking.AMPConfig,
    pb.endgame.PowerSeriesConfig,
])
def test_config_struct_round_trips(cls):
    cfg = cls()
    r = pickle.loads(pickle.dumps(cfg))
    assert r == cfg
    assert copy.deepcopy(cfg) == cfg


# ----------------------------- result metadata -----------------------------

@pytest.mark.parametrize("cls", [
    pb.nag_algorithm.SolutionMetaDataDoublePrec,
    pb.nag_algorithm.SolutionMetaDataMultiPrec,
])
def test_solution_metadata_round_trips_including_enum_field(cls):
    m = cls()
    m.endgame_success_code = pb.SuccessCode.GoingToInfinity
    m.condition_number = 12.5
    r = pickle.loads(pickle.dumps(m))
    assert r.endgame_success_code == m.endgame_success_code
    assert r == m
