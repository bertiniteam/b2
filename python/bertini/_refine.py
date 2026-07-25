"""Standalone Newton refinement -- the sharpening primitive.

Wraps the core's tracker-free ``NewtonRefine``: iterate Newton's method on a square
system until consecutive approximations agree to a requested tolerance, at the
precision of the inputs.  The intended customers are DEFLATED systems at singular
points, where deflation restores the quadratic convergence the plain system loses.
"""

import numpy as np

import _pybertini

from .multiprec import complex_mp


def newton_refine(system, start, tolerance=1e-30, max_iterations=30, time=None):
    """Newton-refine a point against a square system, without a tracker.

    Iterates full Newton steps until the infinity norm of a step falls at or below
    ``tolerance`` or ``max_iterations`` steps have been taken.  Evaluation runs at
    the precision of the supplied system and point -- lift both (and raise
    ``bertini.default_precision``) to sharpen beyond their current precision.

    Refining a SINGULAR point requires the system to be its deflation (else Newton
    converges only linearly and the requested tolerance is out of reach); square an
    overdetermined deflated system by randomization first, exactly as the tracking
    layer does.

    Parameters
    ----------
    system : bertini.System
        The square system to refine against (total functions, including patches,
        equal to variables).
    start : array_like
        The approximate solution, length equal to the system's variable count.
    tolerance : float, optional
        Stop when the infinity norm of a Newton step is at or below this.
    max_iterations : int, optional
        Refuse to iterate more than this many times.
    time : complex, optional
        Path-variable value for non-autonomous systems; must be omitted (or None)
        for autonomous systems, which is the ordinary sharpening case.

    Returns
    -------
    point : numpy.ndarray
        The refined point (the best iterate reached, even on failure).
    code : bertini.tracking.SuccessCode
        ``Success`` when the tolerance was reached; ``FailedToConverge`` or
        ``MatrixSolveFailure`` otherwise.
    achieved : float
        Infinity norm of the last Newton step -- the agreement actually achieved.
    iterations : int
        Number of Newton iterations taken.
    """
    arr = np.asarray([c if isinstance(c, complex_mp) else complex_mp(c)
                      for c in np.asarray(start).ravel()])
    if time is None:
        point, code, achieved, iterations = _pybertini.newton_refine(
            system, arr, float(tolerance), int(max_iterations))
    else:
        point, code, achieved, iterations = _pybertini.newton_refine(
            system, arr, float(tolerance), int(max_iterations),
            complex_mp(time))
    return np.asarray(point), code, achieved, iterations
