"""bertini's list containers (ListOfInt, VariableGroup, the solution lists, ...) compare by
value with any Python sequence, so `sys.degrees() == [2, 3]` reads the way it should."""

import pytest

import bertini


@pytest.fixture
def xy():
    return bertini.Variable('x'), bertini.Variable('y')


@pytest.fixture
def system(xy):
    x, y = xy
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x**2 - 3)
    sys.add_function(y**3 - x)
    return sys


def test_degrees_equal_a_list_of_ints(system):
    assert system.degrees() == [2, 3]
    assert [2, 3] == system.degrees()          # the list on the left defers to the container
    assert system.degrees() == (2, 3)
    assert not (system.degrees() != [2, 3])


def test_degrees_differ_by_value_or_length(system):
    assert system.degrees() != [3, 2]
    assert system.degrees() != [2]
    assert system.degrees() != [2, 3, 4]
    assert not (system.degrees() == [3, 2])


def test_two_containers_with_equal_contents_are_equal(system):
    # two calls return two objects; they were unequal when == meant identity
    assert system.degrees() == system.degrees()


def test_non_sequences_are_never_equal(system):
    degrees = system.degrees()
    assert degrees != 2
    assert degrees != None                     # noqa: E711 -- testing == itself
    assert degrees != "23"                     # a string is a sequence, but not of ints
    assert degrees != {2, 3}                   # a set is not a sequence


def test_variable_group_equals_its_variables(xy):
    x, y = xy
    group = bertini.VariableGroup([x, y])
    assert group == [x, y]
    assert group != [y, x]


def test_solution_lists_compare_their_vectors_whole(system):
    solver = bertini.ZeroDimSolver(system)
    solver.solve()
    solutions = solver.all_solutions()
    # each element is a numpy vector, whose == is elementwise; the container asks for all of it
    assert solutions == [s for s in solutions]
    assert solutions != [s for s in solutions][:-1]


def test_containers_are_not_hashable(system):
    # compared by value and mutable, like Python's own list
    with pytest.raises(TypeError):
        hash(system.degrees())
