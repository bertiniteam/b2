"""Getting and setting the random state (bertiniteam/b2#518): the draws after set_random_state are
the draws that followed get_random_state, so a computation written out between its stages continues
with the numbers an uninterrupted run would have drawn."""
import pickle

import pytest

import bertini
from bertini.random import RandomState, get_random_state, set_random_state


def _draws_from_here():
    """A homotopy's worth of draws of every kind: vectors, a solve seed, a matrix."""
    return (repr(list(bertini.random.random_vector(3))),
            bertini.random.derive_solve_seed(),
            repr(bertini.random.random_matrix(2, 3, orthonormal=False).tolist()),
            repr(bertini.random.random_real()))


def test_a_state_set_again_draws_what_followed_getting_it():
    bertini.random.set_random_seed(1042)
    _draws_from_here()
    state = get_random_state()
    after_getting_it = _draws_from_here()

    bertini.random.set_random_seed(77)
    _draws_from_here()
    set_random_state(state)
    assert _draws_from_here() == after_getting_it
    assert bertini.random.get_random_seed() == 1042


def test_a_solve_after_setting_a_state_finds_what_the_solve_after_getting_it_found():
    x, y = bertini.variables(['x', 'y'])
    system = bertini.System()
    system.add_variable_group(x, y)
    system.add_functions([x**2 + y**2 - 1, x - y])

    def solved():
        solver = bertini.nag_algorithm.ZeroDimSolver(system)
        solver.solve()
        return sorted(repr([complex(c) for c in s]) for s in solver.all_solutions())

    bertini.random.set_random_seed(3)
    _draws_from_here()
    state = get_random_state()
    after_getting_it = (solved(), _draws_from_here())

    _draws_from_here()
    set_random_state(state)
    assert (solved(), _draws_from_here()) == after_getting_it


def test_a_state_holds_the_seed_and_ordinal_and_compares_by_its_fields():
    bertini.random.set_random_seed(5)
    bertini.random.derive_solve_seed()
    state = get_random_state()
    assert (state.global_seed, state.seed_chosen, state.solve_ordinal) == (5, True, 1)
    assert get_random_state() == state
    _draws_from_here()
    assert get_random_state() != state


def test_a_state_prints_its_seed_and_ordinal():
    bertini.random.set_random_seed(5)
    bertini.random.derive_solve_seed()
    assert repr(get_random_state()) == '<RandomState: global seed 5, solve ordinal 1>'
    assert repr(RandomState()) == '<RandomState: no seed chosen>'


def test_a_state_built_field_by_field_draws_as_one_got():
    bertini.random.set_random_seed(11)
    _draws_from_here()
    got = get_random_state()
    after_getting_it = _draws_from_here()

    built = RandomState()
    built.global_seed = got.global_seed
    built.seed_chosen = got.seed_chosen
    built.solve_ordinal = got.solve_ordinal
    built.draw_stream = got.draw_stream
    bertini.random.set_random_seed(12)
    set_random_state(built)
    assert _draws_from_here() == after_getting_it


def test_a_state_goes_through_text_and_pickle():
    bertini.random.set_random_seed(13)
    _draws_from_here()
    state = get_random_state()
    text = state.to_text()
    assert text.startswith('b2randstate/1\nglobal_seed 13\nseed_chosen 1\nsolve_ordinal 1\n')
    assert RandomState.from_text(text) == state
    assert pickle.loads(pickle.dumps(state)) == state


def test_a_text_that_is_not_a_state_is_refused():
    bertini.random.set_random_seed(9)
    text = get_random_state().to_text()
    for bad in ('', 'global_seed 9', text.replace('b2randstate/1', 'b2randstate/2'),
                text.split('draw_stream')[0]):
        with pytest.raises(ValueError, match='b2randstate/1'):
            RandomState.from_text(bad)
