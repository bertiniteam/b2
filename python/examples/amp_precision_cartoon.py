"""Recreate the homotopy-continuation "cartoon" from real tracked data, and use it to *see*
adaptive-multiprecision (AMP) precision changes.

The textbook picture of homotopy continuation (see
``doc_resources/images/homotopycontinuation_generic.png``) draws paths flowing from smooth start
points at ``t=1`` (right) to the solutions of the target at ``t=0`` (left), past an "endgame
boundary".  Those drawings are cartoons.  This script makes the *same* picture from a real solve:
we attach a ``PathCollectionObserver`` to the tracker, solve cyclic-5 with the adaptive-precision
Cauchy solver, collect every path into a pandas DataFrame, and plot it.

What the plot shows, from real data:
  * x-axis: ``log10|t|``, so ``t=1`` is on the right and ``t -> 0`` is on the left (log time);
  * height: the real part of a dehomogenized coordinate (paths weave and fan out to the endpoints;
    paths going to infinity shoot off — the homogenizing coordinate goes to zero);
  * color: the condition number along the path (log scale);
  * a vertical line at the endgame boundary (``t = 0.1``);
  * a marker wherever the tracker *raised its working precision* — the thing this whole
    investigation was about.  After the Criterion-B fix (ADR-0038) almost every cyclic-5 path
    stays in double; only a couple, which track unusually deep into the endgame, ever escalate.

Run (needs matplotlib + pandas):
    python python/examples/amp_precision_cartoon.py            # writes SVG + PNG next to cwd
    python python/examples/amp_precision_cartoon.py /tmp/out   # writes into /tmp/out
"""
import os
import sys

import numpy as np
import bertini as pb
from bertini.tracking import observers


N = 5
ENDGAME_BOUNDARY = 0.1   # default t at which the endgame begins
HOMVAR_INDEX = 0         # homogenize() prepends the homogenizing coordinate


def cyclic_system(n):
    """The cyclic-n system (same construction as examples/solve_cyclic.py)."""
    x = [pb.Variable('x{}'.format(i)) for i in range(n)]
    w = x + x
    sys = pb.System()
    for length in range(1, n):
        sys.add_function(np.sum([np.prod(w[start:start + length]) for start in range(n)]))
    sys.add_function(np.prod(x) - 1)
    sys.add_variable_group(pb.VariableGroup(x))
    return sys


def collect_paths():
    """Solve cyclic-5 with AMP+Cauchy, capturing every solution path via an observer.

    Uses the purpose-built SolutionPathCollector: attached to the SOLVER, it grabs the tracker that
    actually runs each path (clone-safe under threading) and collects the whole journey to t -> 0,
    endgame included -- one PathDataCollector per solution path.  Returns a list of per-path
    DataFrames.
    """
    pb.random.set_random_seed(1)   # deterministic homotopy
    solver = pb.nag_algorithm.ZeroDimCauchyAdaptivePrecisionTotalDegree(cyclic_system(N))

    collector = pb.nag_algorithm.observers.SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()

    return [s.as_dataframe() for s in collector.series if len(s) > 0]


def make_plot(paths, out_stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm
    import matplotlib.cm as cm

    fig, ax = plt.subplots(figsize=(11, 6.5))

    # global condition-number range, for a shared color scale
    all_cond = np.concatenate([p["condition_number"].to_numpy() for p in paths])
    all_cond = all_cond[all_cond > 0]
    norm = LogNorm(vmin=max(all_cond.min(), 1.0), vmax=all_cond.max())
    cmap = matplotlib.colormaps["viridis"]

    esc_x, esc_y = [], []   # precision-increase markers
    for p in paths:
        t = p["t"].to_numpy()
        abst = np.abs(t)
        good = abst > 0
        x = np.log10(abst[good])
        z = p["z{}".format(HOMVAR_INDEX + 1)].to_numpy()[good]
        hom = p["z{}".format(HOMVAR_INDEX)].to_numpy()[good]
        # dehomogenize; height = real part of an affine coordinate
        with np.errstate(divide="ignore", invalid="ignore"):
            y = np.real(z / hom)
        cond = p["condition_number"].to_numpy()[good]
        prec = p["precision"].to_numpy()[good]

        finite = np.isfinite(x) & np.isfinite(y)
        x, y, cond, prec = x[finite], y[finite], cond[finite], prec[finite]
        if len(x) < 2:
            continue

        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=1.0, alpha=0.75)
        lc.set_array(np.clip(cond[:-1], norm.vmin, norm.vmax))
        ax.add_collection(lc)

        inc = np.where(np.diff(prec) > 0)[0] + 1   # indices where precision rose
        esc_x.extend(x[inc]); esc_y.extend(y[inc])

    if esc_x:
        ax.scatter(esc_x, esc_y, marker="*", s=55, facecolor="crimson",
                   edgecolor="black", linewidths=0.3, zorder=5, alpha=0.85,
                   label="precision raised (-> multiprecision)")

    ax.axvline(np.log10(ENDGAME_BOUNDARY), color="purple", ls="--", lw=1.5)
    ax.text(np.log10(ENDGAME_BOUNDARY), ax.get_ylim()[1], "  endgame boundary",
            color="purple", va="top", ha="left", fontsize=9)

    ax.set_xlabel(r"$\log_{10}|t|$   (start $t=1$ at right $\;\to\;$ target $t=0$ at left)")
    ax.set_ylabel(r"$\mathrm{Re}(x_1)$  (dehomogenized)")
    ax.set_title("cyclic-5 homotopy paths from real data — colored by condition number\n"
                 "(adaptive precision; stars mark where the tracker left double precision)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.margins(x=0.02)
    if esc_x:
        ax.legend(loc="lower left", fontsize=9)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax); cb.set_label("condition number")

    fig.tight_layout()
    svg, png = out_stem + ".svg", out_stem + ".png"
    fig.savefig(svg);  fig.savefig(png, dpi=150)
    plt.close(fig)
    return svg, png


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    paths = collect_paths()

    n_escalating = sum(int((np.diff(p["precision"].to_numpy()) > 0).any()) for p in paths)
    print("collected {} paths; {} ever raised precision".format(len(paths), n_escalating))

    stem = os.path.join(out_dir, "amp_precision_cartoon_cyclic5")
    svg, png = make_plot(paths, stem)
    print("wrote {}\n      {}".format(svg, png))


if __name__ == "__main__":
    main()
