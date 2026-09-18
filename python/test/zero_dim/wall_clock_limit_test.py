# This file is part of Bertini 2.
#
# python/test/zero_dim/wall_clock_limit_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/zero_dim/wall_clock_limit_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/zero_dim/wall_clock_limit_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""A wall-clock budget on paths, from Python: the solver knob, the bare tracker's deadline, and the
recall policy's reading of an abandonment.  Correctness lives in C++; this is the interface.

Budgets in these tests are either a nanosecond (gone before the first step, on any machine) or an
hour (never reached), so nothing here depends on how fast the machine is.
"""

import bertini as pb
from bertini import ZeroDimSolver
from bertini.nag_algorithm import ZeroDimConfig, RecordsConfig, RecallPolicy


def _two_quadrics():
    x, y = pb.variables(['x', 'y'])
    sys = pb.System()
    sys.add_variable_group(x, y)
    sys.add_functions([x * x - 1, y * y - 1])
    return sys


def test_a_budget_no_path_can_meet_abandons_every_path_with_its_stamp():
    solver = ZeroDimSolver(_two_quadrics())
    solver.update(max_path_wall_clock_duration=1e-9)
    solver.solve()

    metas = solver.solution_metadata()
    assert len(metas) == 4
    for m in metas:
        assert m.pre_endgame_success_code == pb.SuccessCode.WallClockLimitReached
        assert len(m.last_point) > 0
        assert m.wall_clock_limit_seconds == 1e-9
    assert not solver.was_stopped_early()          # ran out of patience, not interrupted
    assert not solver.get_tracker().has_max_wall_clock_time()   # the deadline was the path's, not left behind


def test_a_generous_budget_changes_nothing():
    solver = ZeroDimSolver(_two_quadrics())
    solver.update(max_path_wall_clock_duration=3600)
    solver.solve()
    for m in solver.solution_metadata():
        assert m.endgame_success_code == pb.SuccessCode.Success
        assert m.wall_clock_limit_seconds == 3600
    assert solver.get_config(ZeroDimConfig).max_path_wall_clock_duration == 3600


def test_a_whole_solve_budget_reads_as_an_interrupt():
    solver = ZeroDimSolver(_two_quadrics())
    solver.update(max_solve_wall_clock_duration=1e-9)
    solver.solve()                                  # returns: a budget is not an error
    assert solver.was_stopped_early()
    assert solver.num_paths_never_started() == 4
    assert all(m.pre_endgame_success_code == pb.SuccessCode.NeverStarted
               for m in solver.solution_metadata())

    solver.update(max_solve_wall_clock_duration=3600)
    solver.solve()
    assert not solver.was_stopped_early()
    assert all(m.endgame_success_code == pb.SuccessCode.Success for m in solver.solution_metadata())


def test_a_bare_tracker_carries_a_deadline():
    tracker = ZeroDimSolver(_two_quadrics()).get_tracker()
    assert not tracker.has_max_wall_clock_time()
    tracker.set_max_wall_clock_duration(3600)
    assert tracker.has_max_wall_clock_time()
    tracker.clear_max_wall_clock_time()
    assert not tracker.has_max_wall_clock_time()


def test_recall_policy_reads_an_abandonment_under_a_budget(tmp_path):
    rec = str(tmp_path / 'records')

    def solve(limit, policy):
        pb.random.set_random_seed(42)
        solver = ZeroDimSolver(_two_quadrics())
        solver.record_to(rec)
        solver.update(max_path_wall_clock_duration=limit, recall=policy)
        solver.solve()
        return solver

    a = solve(1e-9, RecallPolicy.Completed)
    assert a.num_paths_recalled() == 0
    b = solve(1e-9, RecallPolicy.Completed)         # same patience: the abandonments stand in
    assert b.num_paths_recalled() == 4
    assert all(m.pre_endgame_success_code == pb.SuccessCode.WallClockLimitReached
               for m in b.solution_metadata())
    c = solve(3600, RecallPolicy.Everything)        # told to reuse regardless
    assert c.num_paths_recalled() == 4
    d = solve(3600, RecallPolicy.Completed)         # more patience: re-tracked, and they finish
    assert d.num_paths_recalled() == 0
    assert all(m.endgame_success_code == pb.SuccessCode.Success for m in d.solution_metadata())
