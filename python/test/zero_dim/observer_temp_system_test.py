"""Regression for ADR-0039: solving with an observer must not crash when the System was passed as
a temporary (e.g. built in a helper function and not kept in a variable).

Before the fix, the ZeroDim binding did not keep the given System alive: the solver's internal
(shallow) copy shared function-tree nodes whose shared_ptr deleters are tied to the Python wrappers,
so once the temporary System was collected, the endgame's Differentiate() dropped a node and the
deleter touched a freed Python object -> SIGSEGV.  The binding now ties the System's lifetime to the
solver (with_custodian_and_ward), so this idiomatic usage is safe.
"""
import bertini as pb
from bertini.nag_algorithm import ZeroDim, observers as nobs


def _build_small_system():
    """A small system, built in a helper so its Variable/System wrappers are temporaries.

    f1 = (x*y - 3x + 2)*(x - 4),  f2 = y - x^2  -> a singular endpoint at (1,1), nonsingular at
    (-2,4) and (4,16), and paths that diverge to infinity (6 total).
    """
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function((x * y - 3 * x + 2) * (x - 4))
    sys.add_function(y - x ** 2)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    return sys


def test_observer_with_temporary_system_does_not_crash():
    pb.random.set_random_seed(1)
    # System built in a function and passed straight in -- no Python reference retained.
    zd = ZeroDim(_build_small_system(), mptype='adaptive')
    collector = nobs.SolutionPathCollector()
    zd.add_observer(collector)
    zd.solve()                                   # must not segfault

    assert len(zd.finite_solutions()) == 4       # (1,1) [singular], (-2,4), (4,16)
    assert len(collector.series) == 6            # every path collected
    # the collected per-path data is usable after the solve (touches the shared nodes again)
    for series in collector.series:
        assert len(series) >= 1
        series.as_dataframe()
