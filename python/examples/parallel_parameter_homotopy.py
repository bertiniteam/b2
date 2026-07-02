#!/usr/bin/env python3
"""Parallel parameter homotopy: a chemical reaction network's bistability map.

Two-level parallelism for a parameter sweep:

* **MPI across parameter points** -- each rank owns a slice of the parameter grid (this is the
  *outer* loop; the points are independent, so no manager-worker dispatch is needed -- just split
  the grid and gather the answers).
* **threads across paths** -- inside each point's solve, the path tracking runs on all the cores
  the rank is given (``OMP_NUM_THREADS``).

The model is the **Schlogl** reaction network, whose steady state is a cubic in the species
concentration ``X``::

    k2 * X^3  -  k1 * X^2  +  k4 * X  -  k3  =  0

Over the positive-orthant parameter space it is **bistable**: depending on the rate constants it
has either **one** or **three** positive real steady states.  We map that region by counting the
positive real solutions at each grid point -- letting the solver do the real/finite classification,
never a hand-rolled cutoff.

Run it (per the MPI notes in the "Solving at scale" tutorial -- always run a *file*, never a
heredoc; ``--bind-to none`` lets each rank's threads spread across cores)::

    python parallel_parameter_homotopy.py                       # serial baseline
    mpirun -n 4 python parallel_parameter_homotopy.py           # 4 ranks, 1 thread each
    OMP_NUM_THREADS=3 mpirun -n 4 --bind-to none \
        python parallel_parameter_homotopy.py                   # 4 ranks x 3 threads

The result (a counts grid, and optionally a heatmap PNG) is assembled on rank 0.
"""

import argparse

import numpy as np

import bertini as pb
from bertini.nag_algorithm import parameter_sweep

# The species variable is SHARED across every member of the family: a parameter homotopy
# interpolates the members' coefficients, so they must be built over the same Variable object.
X = pb.Variable('X')


def make_system(params):
    """The Schlogl steady-state system at the given rate constants (k1, k2, k3, k4)."""
    k1, k2, k3, k4 = params
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([X]))
    sys.add_function(k2 * X * X * X - k1 * X * X + k4 * X - k3)
    return sys


def real_coeff(value):
    """An exact real coefficient node (a complex literal with zero imaginary part)."""
    return pb.multiprec.complex_mp(repr(float(value)), "0")


def positive_real_count(solver):
    """Number of positive real steady states, from the solver's OWN classification.

    ``to_dataframe`` tags each endpoint is_real / is_finite (and multiplicity, etc.); we keep the
    real, finite ones and require the concentration to be positive -- a CRN steady state is a
    positive real point.  We do not re-derive real-ness ourselves.
    """
    df = solver.to_dataframe()
    n = 0
    for _, row in df.iterrows():
        if not (row['is_real'] and row['is_finite']):
            continue
        sol = row['solution']
        x = complex(sol[0] if hasattr(sol, '__len__') else sol)
        if x.real > 1e-9:
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grid', type=int, default=24, help='grid resolution per axis (default 24)')
    ap.add_argument('--save', default=None, help='save the heatmap to this PNG/SVG path (rank 0)')
    args = ap.parse_args()

    # MPI is optional: with mpi4py absent (or one process) this is a plain serial/threaded sweep.
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank, size = comm.Get_rank(), comm.Get_size()
    except ImportError:
        comm, rank, size = None, 0, 1

    N = args.grid
    k1s = np.linspace(3.5, 8.5, N)        # the X^2 rate
    k3s = np.linspace(1.0, 13.0, N)       # the influx term
    # one generic (complex) member, solved once; its solutions seed every track
    generic = (pb.multiprec.complex_mp("0.7", "0.4"),  pb.multiprec.complex_mp("1", "0.3"),
               pb.multiprec.complex_mp("0.31", "-0.52"), pb.multiprec.complex_mp("1.2", "0.9"))
    # the parameter grid, flattened row-major into one list of targets
    targets = [(real_coeff(k1), real_coeff(1), real_coeff(k3), real_coeff(11))
               for k3 in k3s for k1 in k1s]

    # The entire two-level parallel sweep -- MPI across the grid points, threads across the paths
    # inside each point -- is this one call.  `comm=` spreads the points over the ranks and gathers
    # the results; `collect=` reduces each solve to a picklable number before it crosses ranks.
    # (Drop `comm=` for a plain serial/threaded run.)
    counts_flat = parameter_sweep(make_system, generic, targets,
                                  collect=positive_real_count, comm=comm)

    if rank == 0:
        counts = np.array(counts_flat, dtype=int).reshape(N, N)
        seen = sorted(set(counts.ravel().tolist()))
        print("swept %d points over %d rank(s); positive-real counts seen: %s"
              % (N * N, size, seen))
        print("bistable (3 positive real) at %d of %d grid points"
              % (int((counts == 3).sum()), N * N))

        if args.save:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 5))
            im = ax.imshow(counts, origin='lower', aspect='auto', cmap='RdYlBu_r',
                           extent=[k1s[0], k1s[-1], k3s[0], k3s[-1]], vmin=1, vmax=3)
            cb = fig.colorbar(im, ax=ax, ticks=[1, 2, 3])
            cb.set_label('positive real steady states')
            ax.set_xlabel(r'$k_1$ (X$^2$ rate)'); ax.set_ylabel(r'$k_3$ (influx)')
            ax.set_title('Schlögl model: bistability region')
            fig.tight_layout(); fig.savefig(args.save)
            print("saved", args.save)


if __name__ == '__main__':
    main()
