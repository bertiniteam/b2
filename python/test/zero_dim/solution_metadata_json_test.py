"""A solution's metadata as exact JSON, and back (bertiniteam/b2#521): every field reads back
equal, so a program that keeps solutions between runs keeps what b2 knew about them."""
import copy
import pickle

import pytest

import bertini as pb
from bertini import ZeroDimSolver

FIELDS = ('path_index', 'solution_index', 'precision_changed', 'time_of_first_prec_increase',
          'max_precision_used', 'path_time_seconds', 'num_successful_steps', 'num_failed_steps',
          'wall_clock_limit_seconds', 'pre_endgame_success_code', 'condition_number',
          'newton_residual', 'final_time_used', 'accuracy_estimate_internal_coords',
          'accuracy_estimate_user_coords', 'cycle_num', 'precision_digits', 'accuracy_digits',
          'endgame_success_code', 'function_residual', 'multiplicity',
          'multiplicity_representative', 'is_real', 'is_finite', 'is_singular',
          'crossing_unresolved', 'is_nonsolution')


def _metadata(mptype):
    pb.random.set_random_seed(4)
    x, y = pb.variables(['x', 'y'])
    system = pb.System()
    system.add_variable_group([x, y])
    system.add_functions([x ** 2 + y ** 2 - 1, x - y])
    solver = ZeroDimSolver(system, mptype=mptype)
    solver.solve(show_progress=False)
    return list(solver.solution_metadata())


def _same(a, b):
    for name in FIELDS:
        assert repr(getattr(a, name)) == repr(getattr(b, name)), name
    for name in ('singular_values', 'latest_path_point'):
        assert repr(list(getattr(a, name))) == repr(list(getattr(b, name))), name


@pytest.mark.parametrize('mptype', ['double', 'adaptive'])
def test_a_solves_metadata_reads_back_equal_in_every_field(mptype):
    for m in _metadata(mptype):
        text = m.to_json()
        metadata_loaded = type(m).from_json(text)
        _same(metadata_loaded, m)
        assert metadata_loaded.to_json() == text


@pytest.mark.parametrize('mptype', ['double', 'adaptive'])
def test_metadata_pickles_and_copies_whole(mptype):
    """Every field, the read-only ones too: the singular values and the point a path stopped at
    were dropped when pickling went through the writable fields."""
    metadata = _metadata(mptype)
    assert any(len(m.singular_values) for m in metadata)
    for m in metadata:
        for metadata_loaded in (pickle.loads(pickle.dumps(m)), copy.copy(m), copy.deepcopy(m)):
            _same(metadata_loaded, m)


def test_text_that_is_not_such_metadata_is_refused():
    m = _metadata('double')[0]
    with pytest.raises(ValueError, match='b2meta/1'):
        type(m).from_json('{"format": "b2meta/2"}')
    with pytest.raises(ValueError):
        type(m).from_json('not json')
