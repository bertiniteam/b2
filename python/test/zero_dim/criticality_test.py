"""Critical points of a space curve via the rank-deficient (nullvector) Jacobian.

Critical points of a curve f=0 with respect to a linear projection pi are where the projection's
differential vanishes on the curve's tangent.  The degree-smart way to state this -- instead of a
high-degree determinant -- is that the Jacobian of f stacked with pi's coefficient row is
**rank-deficient**: the 3x3 matrix M = [ J_f ; pi ] has a nonzero null vector v.  We solve for
(x, y, z) on the curve together with that null vector, plus one patch equation so v != 0.

The curve is reducible (Bezout degree 3*3 = 9): a line (the y-axis), two interlocking circles, and
a degree-4 quartic where the two cylinders meet.  The criticality system finds critical points on
all of it -- smooth projection-critical points plus the singular crossings where components meet.
This test checks the *checkable subset*: each circle is genus-0 degree-2, so its two critical
points are predictable in closed form -- the analytic oracle that tells us the solve is right.

    f = x (x^2 + y^2 - 1)          # plane x=0  u  cylinder x^2 + y^2 = 1
    g = z ((y-1)^2 + z^2 - 1)      # plane z=0  u  cylinder (y-1)^2 + z^2 = 1
"""

import numpy as np
import pytest

import bertini as pb


def _curve():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    f = x * (x**2 + y**2 - 1)
    g = z * ((y - 1)**2 + z**2 - 1)
    return x, y, z, f, g


def _criticality_system(x, y, z, f, g, pi):
    """The 6x6 nullvector criticality system for the projection with coefficient row ``pi``."""
    J = pb.jacobian([f, g], [x, y, z])                       # 2 x 3 symbolic
    M = np.vstack([J, pb.coefficients(pi)])           # 3 x 3: J_f over the projection gradient
    v = np.array(pb.variables('v', 3), dtype=object)                       # the null-vector unknowns

    sys = pb.System()
    sys.add(pb.VariableGroup([x, y, z, *v]), f, g)           # curve equations
    sys.add_functions(M @ v)                         # M v = 0   (3 equations)
    patch = pb.random_matrix(1, 3, symbolic=True)            # h . v = 1 keeps v away from zero
    sys.add_function((patch @ v)[0] - 1)
    return sys


def test_critical_points_of_interlocking_circles():
    pb.random.set_random_seed(165)
    x, y, z, f, g = _curve()

    # a real projection, so the two real circles' critical points are real and in closed form
    pi = pb.random_matrix(1, 3, real=True, orthonormal=False)
    a, b, c = (complex(e).real for e in pi.ravel())          # pi = (a, b, c)

    sys = _criticality_system(x, y, z, f, g, pi)
    assert list(sys.degrees()) == [3, 3, 3, 3, 1, 1]         # the criticality system, 6 eqns / 6 vars

    solver = pb.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive',
                                            startsystem='binomial')
    solver.solve()

    finite = []
    for s in solver.solutions():
        p = np.array(s)
        xyz = (complex(p[0]), complex(p[1]), complex(p[2]))
        if all(abs(w) < 1e6 for w in xyz):
            finite.append(xyz)
    assert finite, "expected finite critical points"

    # --- analytic oracle: the two critical points on each real circle ---------------------------
    # circle A in the plane x=0:  (y-1)^2 + z^2 = 1.  pi restricted = b*y + c*z (a*x = 0), extremized
    # at (0, 1 +/- b/r, +/- c/r) with r = sqrt(b^2 + c^2).
    r = np.hypot(b, c)
    # circle B in the plane z=0:  x^2 + y^2 = 1.  pi restricted = a*x + b*y, extremized at
    # (+/- a/s, +/- b/s, 0) with s = sqrt(a^2 + b^2).
    s = np.hypot(a, b)
    predicted = {
        'A1': (0.0, 1 + b / r, c / r),
        'A2': (0.0, 1 - b / r, -c / r),
        'B1': (a / s, b / s, 0.0),
        'B2': (-a / s, -b / s, 0.0),
    }

    def recovered(point):
        return any(max(abs(w.real - q) + abs(w.imag) for w, q in zip(sol, point)) < 1e-5
                   for sol in finite)

    for name, point in predicted.items():
        assert recovered(point), f"missing predicted critical point {name} = {point}"
        # it lies on the curve
        px, py, pz = point
        assert abs(px * (px**2 + py**2 - 1)) < 1e-9
        assert abs(pz * ((py - 1)**2 + pz**2 - 1)) < 1e-9

    # --- every predicted point is genuinely critical: [ J_f(point) ; pi ] is rank-deficient -----
    fg = pb.System()
    fg.add_variable_group(pb.VariableGroup([x, y, z]))
    fg.add_function(f)
    fg.add_function(g)
    for point in predicted.values():
        raw = np.array(fg.eval_jacobian(np.array(point, dtype=complex))).reshape(2, 3)
        Jf = np.array([[complex(e) for e in row] for row in raw])    # to plain complex
        M = np.vstack([Jf, np.array([a, b, c], dtype=complex)])
        assert np.linalg.svd(M, compute_uv=False)[-1] < 1e-9   # rank deficient -> critical here
