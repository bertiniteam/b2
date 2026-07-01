"""Move a slice, hold the system fixed -- Bertini 2 tutorial (moving_slice).

Run:  python moving_slice.py
"""

import numpy as np
import bertini
from bertini import linalg, nag_algorithm


def move_one_slice(gamma):
    """A fixed unit circle sliced by a line sliding from y=0 to y-x=0."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    fixed = bertini.System()
    fixed.add_variable_group(bertini.VariableGroup([x, y]))
    fixed.add_function(x*x + y*y - 1)                 # the fixed unit circle

    start_moving = bertini.System()                   # the slice at t=1: the x-axis y = 0
    start_moving.add_variable_group(bertini.VariableGroup([x, y]))
    start_moving.add_function(y)

    end_moving = bertini.System()                     # the slice at t=0: the diagonal y - x = 0
    end_moving.add_variable_group(bertini.VariableGroup([x, y]))
    end_moving.add_function(y - x)

    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)

    target = bertini.system.concatenate(fixed, end_moving)   # the t=0 system: circle + diagonal
    start_points = [np.array([bertini.multiprec.Complex('1'),  bertini.multiprec.Complex('0')]),
                    np.array([bertini.multiprec.Complex('-1'), bertini.multiprec.Complex('0')])]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4))
                   for s in solver.all_solutions())

    r = round(1 / np.sqrt(2), 4)
    assert roots == sorted([(r, r), (-r, -r)])        # circle ∩ diagonal = (±1/√2, ±1/√2)


def static_and_moving_slice(gamma):
    """A fixed unit sphere cut by a static slice z=0 and a moving slice y=0 -> y-x=0."""
    x, y, z = bertini.Variable('x'), bertini.Variable('y'), bertini.Variable('z')

    fixed = bertini.System()
    fixed.add_variable_group(bertini.VariableGroup([x, y, z]))
    fixed.add_function(x*x + y*y + z*z - 1)                                # the sphere
    linalg.add_linear(fixed, np.array([[0, 0, 1]]), np.array([x, y, z]))   # static slice z = 0

    start_moving = bertini.System(); start_moving.add_variable_group(bertini.VariableGroup([x, y, z]))
    start_moving.add_function(y)                                           # moving slice at t=1
    end_moving = bertini.System(); end_moving.add_variable_group(bertini.VariableGroup([x, y, z]))
    end_moving.add_function(y - x)                                         # moving slice at t=0

    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)
    target = bertini.system.concatenate(fixed, end_moving)
    start_points = [np.array([bertini.multiprec.Complex(str(a)), bertini.multiprec.Complex('0'),
                              bertini.multiprec.Complex('0')]) for a in (1, -1)]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4), round(complex(s[2]).real, 4))
                   for s in solver.all_solutions())
    r = round(1 / np.sqrt(2), 4)
    assert roots == sorted([(r, r, 0.0), (-r, -r, 0.0)])

    # The fixed system is left out of the motion: dH/dt is zero on the fixed blocks.
    pt = np.array([bertini.multiprec.Complex('0.3'),
                   bertini.multiprec.Complex('0.4'),
                   bertini.multiprec.Complex('0.5')])
    # the adaptive solve above left H at double precision; match it to the evaluation point's
    # precision before evaluating the homotopy directly.
    H.precision(pt[0].precision)
    dHdt = H.eval_time_derivative(pt, bertini.multiprec.Complex('0.5'))
    assert abs(complex(dHdt[0])) == 0.0      # sphere row: out of dH/dt
    assert abs(complex(dHdt[1])) == 0.0      # static slice row: out of dH/dt
    assert abs(complex(dHdt[2])) > 0.0       # only the moving slice carries t


def deform_product_into_polynomial(gamma):
    """Deform a product of linears (x-1)(x+1) into the circle, slice y=1/2 static."""
    x, y = bertini.Variable('x'), bertini.Variable('y')

    fixed = bertini.System(); fixed.add_variable_group(bertini.VariableGroup([x, y]))
    linalg.add_linear(fixed, np.array([[0, 1]]), np.array([x, y]), ['-1/2'])    # static slice y = 1/2

    start_moving = bertini.System(); start_moving.add_variable_group(bertini.VariableGroup([x, y]))
    linalg.add_products_of_linears(start_moving, [[[1, 0, -1], [1, 0, 1]]])      # (x-1)(x+1), a structured block

    end_moving = bertini.System(); end_moving.add_variable_group(bertini.VariableGroup([x, y]))
    end_moving.add_function(x*x + y*y - 1)                                       # the polynomial

    H = nag_algorithm.moving_homotopy(fixed, start_moving, end_moving, gamma=gamma)
    target = bertini.system.concatenate(fixed, end_moving)
    start_points = [np.array([bertini.multiprec.Complex(str(a)), bertini.multiprec.Complex('0.5')])
                    for a in (1, -1)]

    solver = nag_algorithm.HomotopySolver(H, start_points, target)
    solver.solve()
    roots = sorted((round(complex(s[0]).real, 4), round(complex(s[1]).real, 4))
                   for s in solver.all_solutions())
    s3 = round(np.sqrt(3) / 2, 4)
    assert roots == sorted([(s3, 0.5), (-s3, 0.5)])     # circle ∩ {y = 1/2}


def main():
    gamma = linalg.coefficient(bertini.multiprec.Complex('0.6', '0.8'))   # off the real axis
    move_one_slice(gamma)
    static_and_moving_slice(gamma)
    deform_product_into_polynomial(gamma)


if __name__ == '__main__':
    main()
