"""Precision models: double, multiple, adaptive -- Bertini 2 tutorial.

The three precision models differ in the interface: the types you get back
and how you set the precision.
Run:  python precision_models.py
"""

import numpy as np
import bertini


def build_system():
    """Build the polynomial system used throughout the tutorial."""
    x, y = bertini.Variable('x'), bertini.Variable('y')
    system = bertini.System()
    system.add_function(x*x + y*y - 1)
    system.add_function(x + y)
    system.add_variable_group(bertini.VariableGroup([x, y]))
    return system


def solver_types(system):
    """mptype selects the solver class (endgame + precision model)."""
    names = {mptype: type(bertini.ZeroDimSolver(system, mptype=mptype)).__name__
             for mptype in ('double', 'multiple', 'adaptive')}
    assert names['double']   == 'ZeroDimSolverCauchyDoublePrecision'
    assert names['multiple'] == 'ZeroDimSolverCauchyFixedMultiplePrecision'
    assert names['adaptive'] == 'ZeroDimSolverCauchyAdaptivePrecision'


def reading_solutions(system):
    """The model changes the type of the numbers; convert with complex()."""
    bertini.random.set_random_seed(2)

    dbl = bertini.ZeroDimSolver(system, mptype='double'); dbl.solve()
    assert dbl.all_solutions()[0].dtype == np.complex128

    amp = bertini.ZeroDimSolver(system, mptype='adaptive'); amp.solve()
    assert str(amp.all_solutions()[0].dtype) == 'complex_mp'     # bertini.complex_mp

    # the same code reads either one:
    def to_python(solution):
        return np.array([complex(c) for c in solution])

    for solver in (dbl, amp):
        pts = sorted(tuple(np.round(to_python(s).real, 4)) for s in solver.all_solutions())
        assert pts == [(-0.7071, 0.7071), (0.7071, -0.7071)]


def setting_precision(system):
    """Fixed multiple works at one precision everywhere; adaptive manages its own."""
    bertini.default_precision(40)                  # 40 digits for this solve
    system.precision(40)                           # the system must match
    m = bertini.ZeroDimSolver(system, mptype='multiple')
    m.solve()
    assert len(m.all_solutions()) == 2

    bertini.default_precision(30)                  # restore a modest default
    system.precision(30)

    from bertini.tracking import AMPConfig
    amp = bertini.ZeroDimSolver(system, mptype='adaptive')
    assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 300
    amp.get_tracker().update(maximum_precision=200)     # tighten the ceiling
    assert amp.get_tracker().get_config(AMPConfig).maximum_precision == 200


def shared_config_surface(system):
    """Apart from precision-specific knobs, the configs are the same across models."""
    shared = {'tolerances', 'zero_dim', 'post_processing', 'auto_retrack'}
    for mptype in ('double', 'multiple', 'adaptive'):
        names = set(bertini.ZeroDimSolver(system, mptype=mptype).config_names())
        assert shared <= names


def main():
    system = build_system()
    solver_types(system)
    reading_solutions(system)
    setting_precision(system)
    shared_config_surface(system)


if __name__ == '__main__':
    main()
