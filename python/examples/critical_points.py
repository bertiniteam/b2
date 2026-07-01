"""Critical points of a space curve via the rank-deficient (nullvector) Jacobian.

The critical points of a curve ``f = 0`` with respect to a linear projection ``pi`` are the points
where ``pi``'s differential vanishes on the curve's tangent.  Rather than a high-degree determinant,
state it as a *rank deficiency*: the Jacobian of ``f`` stacked with ``pi``'s coefficient row is a
3x3 matrix ``M = [ J_f ; pi ]`` with a nonzero null vector ``v``.  Solve for the point ``(x,y,z)``
on the curve together with that null vector, plus one patch equation so ``v != 0``.

The curve here is two interlocking circles (a reducible space curve):

    f = x (x^2 + y^2 - 1)       # plane x=0  u  cylinder x^2 + y^2 = 1
    g = z ((y-1)^2 + z^2 - 1)   # plane z=0  u  cylinder (y-1)^2 + z^2 = 1

Run:  python critical_points.py
"""

import numpy as np

import bertini as pb
from bertini import linalg
from bertini.nag_algorithm import ZeroDimSolver


def main():
    pb.random.set_random_seed(165)
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    f = x * (x**2 + y**2 - 1)
    g = z * ((y - 1)**2 + z**2 - 1)

    # a random real projection; its gradient row is just its (constant) coefficients
    pi = pb.random_matrix(1, 3, real=True, orthonormal=False)

    # --- the three lines that state criticality -------------------------------------------------
    J = pb.jacobian([f, g], [x, y, z])                   # 2 x 3 symbolic Jacobian of the curve
    M = np.vstack([J, linalg.as_coefficients(pi)])       # 3 x 3: J_f stacked over the projection
    v = linalg.variable_vector('v', 3)                   # the null-vector unknowns v0, v1, v2

    sys = pb.System()
    sys.add(pb.VariableGroup([x, y, z, *v]), f, g)       # curve equations
    linalg.add_functions(sys, M @ v)                     # M v = 0   (rank deficiency)
    sys.add_function((pb.random_matrix(1, 3, symbolic=True) @ v)[0] - 1)   # de-zero patch h.v = 1

    solver = ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='binomial')
    solver.solve()

    # keep finite (x, y, z) parts
    points = []
    for s in solver.solutions():
        p = np.array(s)
        xyz = (complex(p[0]), complex(p[1]), complex(p[2]))
        if all(abs(w) < 1e6 for w in xyz):
            points.append(xyz)

    reals = [tuple(round(w.real, 4) for w in p) for p in points if all(abs(w.imag) < 1e-8 for w in p)]
    print(f"{len(points)} finite critical points; {len(reals)} of them real:")
    for p in sorted(set(reals)):
        print("   ", p)


if __name__ == '__main__':
    main()
