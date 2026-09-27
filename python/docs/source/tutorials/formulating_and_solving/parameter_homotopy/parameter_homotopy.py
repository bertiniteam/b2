"""Parameter homotopy: solve once, re-solve many times (Bertini 2 tutorial).

Fixed unit circle meeting a moving horizontal line y = s/2.
Run:  python parameter_homotopy.py
"""

import bertini
from bertini import nag_algorithm

# one line: every solver below records into this directory (see the automatic record keeping
# tutorial).  Pinning the seed makes reruns REPLAY the same homotopies, so a killed
# sweep resumes from the records instead of recomputing.
bertini.records_dir("circle_sweep_records")
bertini.random.set_random_seed(42)

# the variables are SHARED across every member of the family: the parameter homotopy
# interpolates the members' equations, so they must be built over the same Variable objects.
x, y = bertini.Variable('x'), bertini.Variable('y')


def sys_instance(s):
    """The system { x^2 + y^2 - 1, 2y - s }: the unit circle meeting the line y = s/2."""
    sys = bertini.System()
    sys.add_variable_group(bertini.VariableGroup([x, y]))
    sys.add_function(x*x + y*y - 1)
    sys.add_function(2*y - s)
    return sys


def step_one():
    """Pick a generic member and solve it ab initio; its solutions are our start points."""
    start_param_val = bertini.random_complex()
    generic = sys_instance(start_param_val)      
    first = bertini.ZeroDimSolver(generic, mptype='adaptive')
    first.solve()
    start_points = first.all_solutions()
    return generic, start_points


# def step_two(generic, start_points):
#     """Move to another member via a parameter homotopy, without solving from scratch."""
#     target = sys_instance(0)                                 # the line y = 0
#     H = nag_algorithm.straight_line_homotopy(target, generic, gamma=1)
#     solver = bertini.HomotopySolver(H, start_points, target)
#     solver.solve()
#     # moved.all_solutions() are now (+/- 1, 0)
#     return solver


def sweep(generic, start_points):
    """Sweep many parameters, reusing the same start points, never solving from scratch."""
    import numpy as np
    results = {'param_vals':[], 'solns':[]}
    for s in np.linspace(-2, 2, 20, dtype=bertini.real_mp):                               # lines y = 0, -1/2, 1/2
        target = sys_instance(s)
        H = nag_algorithm.straight_line_homotopy(target, generic, gamma=1)
        solver = bertini.HomotopySolver(H, start_points, target)
        solver.solve()
        roots = [p for p in solver.all_solutions()]
        for p in roots:
            xv, yv = complex(p[0]), complex(p[1])
            assert abs(xv*xv + yv*yv - 1) < 1e-8       # on the circle
            assert abs(2*yv - np.float64(s)) < 1e-8                # on the line y = s/2

        results['param_vals'].append(s)
        results['solns'].extend(solver.real_solutions())

    # some converting to facilitate plotting
    results['param_vals'] = np.array(results['param_vals'], dtype=np.float64)
    results['solns'] = np.array(bertini.real(results['solns']), dtype=np.float64)
    return results


def plot(results):
    solns = results['solns']

    import matplotlib.pyplot as plt
    ax = plt.gca()

    line_x_min = -10
    line_x_max = 10

    for s in results['param_vals']:
        # 0 = 2*y - s
        # y = s/2 # a horizontal line at height s/2
        ax.plot([line_x_min, line_x_max], [s/2,s/2],color='xkcd:pale blue',zorder=1) # https://xkcd.com/color/rgb/

    circle = plt.Circle((0,0), 1, color='b', fill=False,zorder=2)
    

    ax.add_patch(circle)
    ax.scatter(solns[:,0],solns[:,1], zorder=3)
    ax.set_aspect('equal', 'box')
    ax.axis([-1.1, 1.1, -1.1, 1.1])

    plt.savefig('parameter_homotopy_circle.png')
    plt.savefig('parameter_homotopy_circle.svg')

def main():
    
    generic, start_points = step_one()
    results = sweep(generic, start_points)

    plot(results)


if __name__ == '__main__':
    main()
