"""Each accuracy estimate says in its name which coordinates it is in (ADR-0069).

``accuracy_estimate_internal_coords`` is the distance between the endgame's last two
approximations in the solver's internal coordinates (homogenized, on the patch): what
``final_tolerance`` is compared with, and a number that does not grow with the solution.
``accuracy_estimate_user_coords`` is the same distance after dehomogenizing: the absolute error in
the units of the caller's variables.  The bare name ``accuracy_estimate`` held the internal one and
read as the user's; it is retired, and says so when read.
"""

import pytest

import bertini as pb
from bertini import ZeroDimSolver

ENDGAMES = ['powerseries', 'cauchy']


def _solve(scale, endgame='powerseries', seed=1):
    """``x^2 = scale^2, y = 2x``: roots at ``(+-scale, +-2 scale)``."""
    pb.random.set_random_seed(seed)
    x, y = pb.variables(['x', 'y'])
    system = pb.System()
    system.add_variable_group([x, y])
    system.add_functions([x ** 2 - pb.coefficient(repr(float(scale * scale))),
                          y - pb.coefficient(2) * x])
    solver = ZeroDimSolver(system, endgame=endgame)
    solver.update(max_norm=1e8)          # the roots are finite; keep the security ceiling above them
    solver.solve(show_progress=False)
    metadata = list(solver.solution_metadata())
    assert [str(m.endgame_success_code) for m in metadata] == ['Success', 'Success']
    return solver, metadata


def test_both_estimates_are_there_under_their_names():
    _, metadata = _solve(3.0)
    for m in metadata:
        assert 0 < float(m.accuracy_estimate_internal_coords) < 1e-8
        assert 0 < float(m.accuracy_estimate_user_coords) < 1e-8


def test_the_bare_name_is_retired_and_says_what_to_read():
    _, metadata = _solve(3.0)
    with pytest.raises(RuntimeError) as refusal:
        metadata[0].accuracy_estimate
    said = str(refusal.value)
    assert 'accuracy_estimate_internal_coords' in said
    assert 'accuracy_estimate_user_coords' in said


def test_the_retirement_is_not_swallowed_by_a_default():
    # getattr with a default swallows an AttributeError and carries on with the default, which
    # is how a reader of the old name would have gone on computing with None
    _, metadata = _solve(3.0)
    with pytest.raises(RuntimeError):
        getattr(metadata[0], 'accuracy_estimate', None)


def test_the_retired_name_is_not_a_field():
    # it is not pickled, printed or tabulated: only the two named estimates are data
    from bertini.config import writable_fields
    _, metadata = _solve(3.0)
    fields = writable_fields(type(metadata[0]))
    assert 'accuracy_estimate_internal_coords' in fields
    assert 'accuracy_estimate_user_coords' in fields
    assert 'accuracy_estimate' not in fields


def test_metadata_round_trips_with_both_estimates():
    import copy
    import pickle
    _, metadata = _solve(3.0)
    m = metadata[0]
    for other in (copy.deepcopy(m), pickle.loads(pickle.dumps(m))):
        assert float(other.accuracy_estimate_internal_coords) == float(m.accuracy_estimate_internal_coords)
        assert float(other.accuracy_estimate_user_coords) == float(m.accuracy_estimate_user_coords)


def test_the_dataframe_has_the_named_columns():
    pytest.importorskip('pandas')
    solver, _ = _solve(3.0)
    columns = list(solver.to_dataframe().columns)
    assert 'accuracy_estimate_internal_coords' in columns
    assert 'accuracy_estimate_user_coords' in columns
    assert 'accuracy_estimate' not in columns


@pytest.mark.parametrize('endgame', ENDGAMES)
def test_the_internal_estimate_does_not_grow_with_the_solution(endgame):
    internal = {scale: max(float(m.accuracy_estimate_internal_coords)
                           for m in _solve(scale, endgame)[1])
                for scale in (1.0, 1e2, 1e4)}
    assert 0.1 < internal[1e4] / internal[1.0] < 10.0, internal
    assert 0.1 < internal[1e4] / internal[1e2] < 10.0, internal


@pytest.mark.parametrize('endgame', ENDGAMES)
def test_the_user_estimate_is_in_the_units_of_the_variables(endgame):
    # ONE problem at two scales under ONE seed: how the two estimates relate to each other
    # depends on the patch that was drawn, so different draws are not comparable
    small = max(float(m.accuracy_estimate_user_coords) for m in _solve(1e2, endgame)[1])
    big = max(float(m.accuracy_estimate_user_coords) for m in _solve(1e4, endgame)[1])
    assert 0.5 < (big / small) / 100.0 < 2.0, (small, big)
