"""Path-crossing detection and the endgame-boundary re-track machinery, at the Python layer.

Two *distinct* paths landing on the same point at the endgame boundary is a probability-0 event:
it is a path crossing (under-resolved tracking), which the zero-dim solver is supposed to detect
and re-track.  Mirrors the C++ tests in ``core/test/nag_algorithms/zero_dim.cpp``; here we exercise
the whole pipeline and the ``endgame_boundary_metadata`` report.

Rather than fish for a crossing inside a big *random* system (slow, and *which* seed crosses depends
on platform floating point), we **construct** a tiny homotopy whose crossing is a deterministic,
portable geometric fact -- a cubic with three solution paths:

    H(x, t) = (x - A(t)) * (x - B) * (x - C),   A(t) = B + delta + kappa*(t - t0)^2,  B=1, C=-2

* ``B`` (flat) and ``C`` (far) are constant paths; ``A(t)`` is a *curved* path that near-misses
  ``B`` from above at ``t0 = 0.9`` (closest gap ``delta``).
* Tracked from ``t=1`` to ``t=0`` with the low-order **Euler** predictor and a *coarse* step, the
  first step lands at ``t=0.9``.  Euler is exact on straight lines but *errs* on a curve, so its
  straight-line prediction of ``A`` overshoots downward, *past* ``B``, into ``B``'s basin.  The
  corrector then lands on ``B`` **exactly** -- a clean root, so the step is *accepted* and adaptive
  step-size control does **not** rescue it -- while ``B`` itself barely moves.

So two paths funnel onto ``B`` and ``A``'s true endpoint ``A(0)`` is lost: a *detectable* crossing.
Finer re-tracking (the solver's resolve step) follows ``A``'s curve correctly and recovers it.  The
overshoot is macroscopic and the coefficients are exact rationals, so the result is identical on
every platform -- no randomness, no seed search, no platform-dependent numerics.  The *fix* is the
solver's re-tracking (and, in real use, a better predictor / tighter tolerance), never tuning; the
detection/resolve logic is also covered by the deterministic C++ tests.
"""

import numpy as np
import pytest

import bertini as pb
from bertini.symbolics import Rational

OK = int(pb.SuccessCode.Success)

_CDT = np.zeros(1, dtype=pb.complex_mp).dtype     # the numpy dtype of a multiprecision complex vector

# Geometry of the planted near-miss (see module docstring).  Exact rationals -> portable.
B_ROOT, C_ROOT = 1, -2                   # the flat path and the far spectator path
T0_NUM, T0_DEN = 9, 10                   # near-miss time t0 = 0.9 = where the first Euler step lands
DELTA_NUM, DELTA_DEN = 1, 10             # closest gap A - B at t0
KAPPA = 10                               # curvature of A: big enough that one Euler step overshoots B
LOST_ROOT = round(B_ROOT + DELTA_NUM / DELTA_DEN + KAPPA * (T0_NUM / T0_DEN) ** 2, 3)  # A(0) = 9.2
KNOWN_DISTINCT = 3                        # A(0), B, C are three distinct solutions at t=0

COARSE_STEP = '0.1'                       # one step from t=1 lands exactly on the near-miss at t=0.9
LOOSE_TOL = 1e-3


def _R(n, d=1):
    return Rational(n, d, 0, 1)


def cubic_crossing_homotopy():
    """Build (homotopy, target, start_points) for the planted cubic crossing."""
    x, t = pb.Variable('x'), pb.Variable('t')
    t0 = _R(T0_NUM, T0_DEN)
    a_of_t = _R(B_ROOT) + _R(DELTA_NUM, DELTA_DEN) + _R(KAPPA) * (t - t0) * (t - t0)

    H = pb.System()
    H.add_variable_group(pb.VariableGroup([x]))
    H.add_function((x - a_of_t) * (x - _R(B_ROOT)) * (x - _R(C_ROOT)))
    H.add_path_variable(t)

    a_at_0 = _R(B_ROOT) + _R(DELTA_NUM, DELTA_DEN) + _R(KAPPA) * (_R(0) - t0) * (_R(0) - t0)
    target = pb.System()
    target.add_variable_group(pb.VariableGroup([x]))
    target.add_function((x - a_at_0) * (x - _R(B_ROOT)) * (x - _R(C_ROOT)))

    a_at_1 = B_ROOT + DELTA_NUM / DELTA_DEN + KAPPA * (1 - T0_NUM / T0_DEN) ** 2
    start_points = [np.array([pb.multiprec.complex_mp(repr(v))], dtype=_CDT)
                    for v in (a_at_1, float(B_ROOT), float(C_ROOT))]
    return H, target, start_points


def _solve(resolve_attempts, predictor=pb.Predictor.Euler, step=COARSE_STEP, tol=LOOSE_TOL):
    """Track the planted cubic and return (report, sorted distinct real roots)."""
    H, target, start_points = cubic_crossing_homotopy()
    solver = pb.nag_algorithm.user_homotopy(H, start_points, target, precision='double')

    solver.get_tracker().predictor(predictor)
    solver.get_tracker().get_stepping().update(initial_step_size=step, max_step_size=step)
    solver.get_tracker().reinitialize_initial_step_size(True)

    tol_cfg = solver.get_config(pb.nag_algorithm.TolerancesConfig)
    tol_cfg.newton_before_endgame = tol
    tol_cfg.newton_during_endgame = tol / 10
    solver.set_config(tol_cfg)

    zd = solver.get_config(pb.nag_algorithm.ZeroDimConfig)
    zd.max_num_crossed_path_resolve_attempts = resolve_attempts
    solver.set_config(zd)

    solver.solve()

    report = solver.endgame_boundary_metadata()
    roots = sorted({round(complex(v[0]).real, 3)
                    for v in solver.all_solutions() if len(v) > 0})
    return report, roots


def test_crossing_is_detected_then_resolved():
    """Report-only loses a solution to the crossing; re-tracking recovers it."""
    # With re-tracking disabled the crossing is detected but left unresolved: two paths funnel onto
    # B=1, so the curved path's true endpoint (A(0)) is lost and the distinct count comes up short.
    report_only, roots_unresolved = _solve(resolve_attempts=0)
    assert not report_only.passed
    assert report_only.num_resolve_attempts == 0
    assert report_only.num_crossings_detected > 0
    assert len(report_only.crossed_path_indices) == report_only.num_crossings_detected
    assert len(roots_unresolved) < KNOWN_DISTINCT
    assert LOST_ROOT not in roots_unresolved          # 9.2 was swallowed by the crossing

    # Re-enable re-tracking: the same crossing is detected, then resolved, and the lost solution is
    # recovered -- the count is now exactly right and A(0) is back.
    resolved, roots_resolved = _solve(resolve_attempts=2)
    assert resolved.num_crossings_detected > 0
    assert resolved.num_resolve_attempts >= 1
    assert resolved.passed
    assert len(roots_resolved) == KNOWN_DISTINCT
    assert roots_resolved == sorted([float(C_ROOT), float(B_ROOT), LOST_ROOT])


def test_clean_solve_reports_no_crossings():
    """With a good (higher-order) predictor the curve is tracked correctly: no crossing at all."""
    report, roots = _solve(resolve_attempts=0, predictor=pb.Predictor.RK4)
    assert report.passed
    assert report.num_crossings_detected == 0
    assert report.num_resolve_attempts == 0
    assert roots == sorted([float(C_ROOT), float(B_ROOT), LOST_ROOT])
