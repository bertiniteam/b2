"""Interface tests for ``bertini.newton_refine`` -- the standalone sharpening
primitive.  Correctness (the isosingular-hierarchy specimens) is gated in C++
(``core/test/nag_algorithms/newton_refine.cpp``); these check the Python surface:
argument handling, return shape, and the friendly refusal of overdetermined
systems."""

import numpy as np
import pytest

import bertini
from bertini.multiprec import complex_mp as C


def _cone_deflated_randomized():
    """The double cone's deflated, hand-randomized (small primes) sharpening system."""
    x, y, z = bertini.variables(list('xyz'))
    S = bertini.System()
    S.add_variable_group([x, y, z])
    f = x**2 + y**2 - z**2
    fx, fy, fz = 2 * x, 2 * y, -2 * z
    S.add_functions([1 * f + 2 * fx + 3 * fy + 5 * fz,
                     7 * f + 11 * fx + 13 * fy + 17 * fz,
                     19 * f + 23 * fx + 29 * fy + 31 * fz])
    return S


@pytest.mark.parametrize("precision", [60], indirect=True)
def test_refines_cone_singularity_on_deflated_system(precision):
    S = _cone_deflated_randomized()
    start = [C('1e-6'), C('-2e-6'), C('5e-7')]
    pt, code, achieved, its = bertini.newton_refine(
        S, start, tolerance=1e-45, max_iterations=50)
    assert str(code) == 'Success'
    assert achieved <= 1e-45
    assert its <= 10                       # quadratic, not limping
    assert isinstance(pt, np.ndarray) and pt.shape == (3,)
    assert max(abs(complex(c)) for c in pt) < 1e-40   # ON the singularity


@pytest.mark.parametrize("precision", [60], indirect=True)
def test_overdetermined_system_raises_with_instructions(precision):
    x = bertini.variables(['x'])[0]
    S = bertini.System()
    S.add_variable_group([x])
    S.add_functions([x**2, x - 1])
    with pytest.raises(RuntimeError, match='[Ss]quare'):
        bertini.newton_refine(S, [C('0.5')], tolerance=1e-20, max_iterations=10)


@pytest.mark.parametrize("precision", [60], indirect=True)
def test_plain_double_root_reports_failure_honestly(precision):
    x = bertini.variables(['x'])[0]
    S = bertini.System()
    S.add_variable_group([x])
    S.add_functions([x**2])
    pt, code, achieved, its = bertini.newton_refine(
        S, [C('1e-6')], tolerance=1e-40, max_iterations=25)
    assert str(code) == 'FailedToConverge'
    assert achieved > 1e-40                # linear halving cannot get there
