"""MPI parallelism support for Bertini2.

Requires an MPI installation (OpenMPI or MPICH) and mpi4py for distributed computing.
When used without MPI, all functions return serial defaults (rank=0, size=1, is_manager=True).
"""

from bertini._pybertini import parallel as _p

rank = _p.rank
size = _p.size
is_manager = _p.is_manager
is_worker = _p.is_worker
initialize = _p.initialize
finalize = _p.finalize

__all__ = ['rank', 'size', 'is_manager', 'is_worker', 'initialize', 'finalize']
