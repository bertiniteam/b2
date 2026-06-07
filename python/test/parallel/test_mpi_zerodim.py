# This file is part of Bertini 2.
#
# python/test/parallel/test_mpi_zerodim.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/parallel/test_mpi_zerodim.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/parallel/test_mpi_zerodim.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""MPI tests for ZeroDim parallel solving.

Tests are skipped automatically if mpi4py is not available.
Run with: mpirun -n <N> python -m pytest python/test/parallel/test_mpi_zerodim.py -v
"""

import pytest

pytest.importorskip("mpi4py")

from mpi4py import MPI
import bertini as pb
from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionTotalDegree


@pytest.fixture
def circle_intersection_solver():
    """x^2 + y^2 - 1 = 0 and x + y = 0: two solutions (±1/√2, ∓1/√2)."""
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_function(x**2 + y**2 - 1)
    sys.add_function(x + y)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    return ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)


def test_parallel_rank_size():
    """Verify bertini.parallel exports rank and size functions."""
    rank = pb.parallel.rank()
    size = pb.parallel.size()
    assert isinstance(rank, int)
    assert isinstance(size, int)
    assert 0 <= rank < size


def test_is_manager_is_worker():
    """Verify is_manager and is_worker are mutually exclusive and cover all ranks."""
    rank = pb.parallel.rank()
    is_mgr = pb.parallel.is_manager()
    is_wrk = pb.parallel.is_worker()

    assert is_mgr ^ is_wrk, f"rank {rank} must be either manager or worker, not both"
    if rank == 0:
        assert is_mgr
    else:
        assert is_wrk


def test_serial_vs_parallel_solution_count(circle_intersection_solver):
    """Solutions in serial match solutions in parallel (rank-0 only).

    All ranks solve the system (required for MPI setup), but only rank 0
    verifies the solution count. Workers have the solutions too but we
    compare counts on the manager.

    For a 2×2 system x^2+y^2-1=0, x+y=0, we expect exactly 2 solutions.
    """
    solver = circle_intersection_solver
    comm = MPI.COMM_WORLD

    solver.solve(communicator=comm)

    if pb.parallel.is_manager():
        solns = solver.solutions()
        assert len(solns) == 2, f"Expected 2 solutions, got {len(solns)}"


def test_serial_solve_works(circle_intersection_solver):
    """Verify serial solve() still works (backward compatibility)."""
    solver = circle_intersection_solver
    solver.solve()

    if pb.parallel.is_manager():
        solns = solver.solutions()
        assert len(solns) == 2
