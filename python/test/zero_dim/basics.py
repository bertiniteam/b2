import pytest

import bertini as pb


@pytest.fixture
def solver():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')

    f1 = x**2 + y**2 - 1
    f2 = x+y

    sys = pb.System()
    sys.add_function(f1)
    sys.add_function(f2)

    sys.add_variable_group(pb.VariableGroup([x, y]))

    return pb.nag_algorithm.ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)


def test_can_solve_multiple_times(solver):
    """
    This test is here because multiple calls to solver.solve() during
    performance testing caused a crash with message

    ```
    Assertion failed: ((bertini::Precision(x(0))==DoublePrecision() || bertini::Precision(x(0)) == Precision()) && "precision of input vector must match current working precision of patch during rescaling"), function RescalePointToFitInPlace, file patch.hpp, line 405.
    ```

    """

    solver.solve()
    solver.solve()
    solver.solve()
