# This file is part of Bertini 2.
#
# python/test/zero_dim/interrupt_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/zero_dim/interrupt_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/zero_dim/interrupt_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Stopping a solve from Python: Ctrl-C, and request_stop() from another thread.

The thing being tested is the seam.  A solve releases the GIL and runs on the calling thread,
so before this a Ctrl-C was noted by CPython and ignored until the solve finished -- the
kernel was trapped and killing the process was the only exit.  Now the solve stops at its
next step, unwinds cleanly, and the solver is left inspectable.

The stop is delivered from INSIDE the solve, by an observer that fires on the n-th completed
path.  A wall-clock timer was tried first and is the wrong tool: a 216-path solve finishes in
under the timer's delay on a fast machine, the signal then lands after the solve's own handler
is gone, and Python's ordinary handler raises -- which looks like a pass and tests nothing.
An observer fires at a point in the computation, not at a time, so it is the same on every
machine.

Two deliveries: `signal.raise_signal(SIGINT)`, which invokes the C-level handler the solve
installed (the Ctrl-C path, which raises KeyboardInterrupt), and `bertini.request_stop()`,
the programmatic path, which does not raise.  raise_signal rather than os.kill because it
behaves identically on every platform the suite runs on.
"""

import signal

import pytest

import bertini as pb


def _many_paths_system(degree=6):
    """A dense three-variable system with degree^3 total-degree paths."""
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y, z]))
    sys.add_function(x**degree + y**degree + z**degree + x*y*z - 1)
    sys.add_function(x**degree - 2*y**degree + z**(degree - 1) + x*y - 3)
    sys.add_function(x**(degree - 1)*z + y**degree - z**degree + x + 2*y - 5)
    return sys


class _StopAfterPaths(pb.nag_algorithm.observers.CustomObserver):
    """Deliver a stop from inside the solve, once `n` paths have completed.

    The observer callback runs with the GIL held, on whichever thread finished the path --
    the calling thread in a serial solve, a worker in a threaded one.  Either is fine for
    both kinds of delivery: the signal handler and request_stop() each do nothing but store
    to an atomic."""

    def __init__(self, n, deliver):
        super().__init__()
        self.n = n
        self.deliver = deliver
        self.completed = 0
        self.fired = False

    def Observe(self, event):
        if isinstance(event, pb.nag_algorithm.observers.PathComplete):
            self.completed += 1
            if self.completed >= self.n and not self.fired:
                self.fired = True
                self.deliver()


def _ctrl_c():
    signal.raise_signal(signal.SIGINT)


def _stopped_taxonomy(solver):
    """Count the kinds of path a stopped solve can leave behind."""
    codes = [m.pre_endgame_success_code for m in solver.solution_metadata()]
    return {
        'never_started': sum(c == pb.SuccessCode.NeverStarted for c in codes),
        'terminated':    sum(c == pb.SuccessCode.ExternallyTerminated for c in codes),
        'success':       sum(c == pb.SuccessCode.Success for c in codes),
        'total':         len(codes),
    }


def test_ctrl_c_during_a_solve_raises_and_leaves_the_solver_inspectable():
    solver = pb.ZeroDimSolver(_many_paths_system(), mptype='adaptive')
    trigger = _StopAfterPaths(3, _ctrl_c)
    solver.add_observer(trigger)
    before = signal.getsignal(signal.SIGINT)

    with pytest.raises(KeyboardInterrupt):
        solver.solve()

    assert trigger.fired                          # the stop really came from inside the solve
    assert solver.was_stopped_early()
    kinds = _stopped_taxonomy(solver)
    assert kinds['never_started'] + kinds['terminated'] > 0
    assert kinds['success'] >= 3                  # at least the paths that triggered it finished
    assert kinds['never_started'] == solver.num_paths_never_started()
    assert kinds['never_started'] + kinds['terminated'] + kinds['success'] == kinds['total']

    # the solver is still a solver: one entry per path, finished ones readable
    assert len(solver.all_solutions()) == kinds['total']

    # Python's own SIGINT handling is exactly as it was, and no request is left standing
    assert signal.getsignal(signal.SIGINT) is before
    assert not pb.stop_requested()


def test_after_an_interrupt_a_plain_ctrl_c_still_reaches_python():
    """The previous handler is restored: once the solve is over, SIGINT is Python's again."""
    solver = pb.ZeroDimSolver(_many_paths_system(), mptype='adaptive')
    solver.add_observer(_StopAfterPaths(3, _ctrl_c))
    with pytest.raises(KeyboardInterrupt):
        solver.solve()

    with pytest.raises(KeyboardInterrupt):
        signal.raise_signal(signal.SIGINT)        # no solve running: plain CPython behaviour


def test_request_stop_from_inside_the_solve_stops_without_raising():
    solver = pb.ZeroDimSolver(_many_paths_system(), mptype='adaptive')
    trigger = _StopAfterPaths(3, pb.request_stop)
    solver.add_observer(trigger)

    result = solver.solve()                       # returns: a programmatic stop is not an error

    assert trigger.fired
    assert solver.was_stopped_early()
    assert not pb.stop_requested()                # withdrawn on the way out
    kinds = _stopped_taxonomy(solver)
    assert kinds['never_started'] + kinds['terminated'] > 0
    assert len(result) <= kinds['success']        # finite solutions come only from finished paths


def test_being_stopped_once_is_not_sticky():
    """The same solver, solved again with nobody stopping it, runs to completion."""
    solver = pb.ZeroDimSolver(_many_paths_system(degree=3), mptype='adaptive')
    trigger = _StopAfterPaths(2, pb.request_stop)
    solver.add_observer(trigger)
    solver.solve()
    assert solver.was_stopped_early()
    solver.remove_observer(trigger)

    solver.solve()
    assert not solver.was_stopped_early()
    assert solver.num_paths_never_started() == 0
    assert all(m.pre_endgame_success_code == pb.SuccessCode.Success
               for m in solver.solution_metadata())


def test_a_stale_request_does_not_stop_the_next_solve():
    """request_stop() means 'stop what is running', not 'do not start': a request made while
    nothing is running is discarded when a solve begins."""
    pb.request_stop()
    assert pb.stop_requested()

    solver = pb.ZeroDimSolver(_many_paths_system(degree=2), mptype='adaptive')
    solver.solve()

    assert not solver.was_stopped_early()
    assert not pb.stop_requested()


def test_homotopy_solver_stops_too():
    """The same seam serves HomotopySolver: it is the same binding and the same tracker check."""
    sys = _many_paths_system(degree=4)
    first = pb.ZeroDimSolver(sys, mptype='adaptive')
    first.solve()
    start = pb.system.start_system.TotalDegreeLinearProduct(sys)
    H = pb.system.make_homotopy(sys, start)
    points = [start.start_point_mp(i) for i in range(start.num_start_points())]

    solver = pb.HomotopySolver(H, points, sys, mptype='adaptive')
    trigger = _StopAfterPaths(2, pb.request_stop)
    solver.add_observer(trigger)
    solver.solve()

    assert trigger.fired
    assert solver.was_stopped_early()
