"""Recreate the *classic homotopy-continuation cartoon* — from real tracked data.

Every textbook draws the same picture (see ``doc_resources/images/homotopycontinuation_generic.png``):
smooth start points at ``t=1`` on the right, paths flowing left to the target's solutions at
``t=0``, past an endgame boundary; the endpoints come in three flavors — nonsingular, singular, and
"at infinity" (paths that diverge).  That drawing is freehand.  This script reproduces it from a real
solve, with each path drawn as a solid line styled by what kind of endpoint it reaches.

We use a small, deliberately-chosen system so there are only a handful of paths and one of each
endpoint flavor:

    f1 = (x*y - 3x + 2)*(x - 4)        f2 = y - x^2

Substituting y = x^2 gives (x-1)^2 (x+2)(x-4) = 0: a DOUBLE root at x=1 (a **singular** endpoint,
reached by two paths) and simple roots at x=-2, 4 (**nonsingular**).  The total-degree start system
has 6 paths, so the remaining two **diverge to infinity** (the homogenizing coordinate -> 0).

Run (needs matplotlib):
    python python/examples/classic_continuation_cartoon.py            # writes SVG + PNG to cwd
    python python/examples/classic_continuation_cartoon.py /tmp/out
"""
import os
import sys

import bertini as pb
from bertini.nag_algorithm import ZeroDim, observers as nobs

ENDGAME_BOUNDARY = 0.1
HOMVAR_INDEX = 0   # homogenize() prepends the homogenizing coordinate; affine coord = z[k]/z[0]
# Height = Re of a fixed GENERIC complex projection of ALL dehomogenized coordinates.  A single
# coordinate's real part collapses the structured total-degree start points (scaled roots of unity)
# onto a few heights; a generic projection (mixing all coordinates with complex weights) separates
# every distinct point.  The diverging paths then run off to large height, so the y-axis is symlog
# (linear near 0 for the finite paths, logarithmic out to the diverging ones).  Fixed for reproducibility.
PROJECTION = (0.6 + 0.8j, -0.9 + 0.4j)

# endpoint-flavor styling (solid lines; color + endpoint marker per flavor), echoing the cartoon
# solid/dashed/dash-dot AND distinct colors per flavor, so the flavors are distinguishable even in
# grayscale / for color-blind readers.
STYLE = {
    "nonsingular": dict(color="#2e7d32", ls="-",  marker="v", label="nonsingular endpoint"),
    "singular":    dict(color="#c2185b", ls="--", marker="*", label="singular endpoint"),
    "infinite":    dict(color="#6a1b9a", ls="-.", marker="^", label="diverges to infinity"),
}


def target_system():
    x, y = pb.Variable("x"), pb.Variable("y")
    sys = pb.System()
    sys.add_function((x * y - 3 * x + 2) * (x - 4))
    sys.add_function(y - x ** 2)
    sys.add_variable_group(pb.VariableGroup([x, y]))
    return sys


def classify(meta):
    if not meta.is_finite:
        return "infinite"
    return "singular" if meta.is_singular else "nonsingular"


def collect():
    """Solve and return [(flavor, DataFrame), ...], one per path."""
    pb.random.set_random_seed(1)
    zd = ZeroDim(target_system(), mptype="adaptive")
    coll = nobs.SolutionPathCollector()
    zd.add_observer(coll)
    zd.solve()
    md = zd.solution_metadata()
    out = []
    for series in coll.series:
        if len(series) == 0:
            continue
        out.append((classify(md[series.path_index]), series.as_dataframe()))
    return out


def make_plot(paths, out_stem):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 6.5))

    def height(df):
        hom = df["z{}".format(HOMVAR_INDEX)].to_numpy()
        n_aff = sum(c.startswith("z") for c in df.columns) - 1   # affine coords (drop homvar)
        with np.errstate(divide="ignore", invalid="ignore"):
            proj = np.zeros(len(df), dtype=complex)
            for k in range(n_aff):
                proj += PROJECTION[k % len(PROJECTION)] * (df["z{}".format(k + 1)].to_numpy() / hom)
        return np.real(proj)                 # generic real projection of the dehomogenized point

    seen = set()
    curves = []     # (flavor, x, y) for every path
    fin_y = []      # heights of the FINITE endpoints, to size the window
    for flavor, df in paths:
        x = np.real(df["t"].to_numpy())
        y = height(df)
        keep = np.isfinite(x) & np.isfinite(y)
        x, y = x[keep], y[keep]
        if len(x) >= 2:
            curves.append((flavor, x, y))
            if flavor != "infinite":
                fin_y.append(y[-1])

    # Window sized to contain the FINITE roots; paths diverging to infinity run off the top, so they
    # do not compress everything else into a flat line (the whole reason for a cutoff).
    lo, hi = (min(fin_y), max(fin_y)) if fin_y else (-1.0, 1.0)
    pad = 0.35 * (hi - lo) + 1.0
    ax.set_ylim(lo - pad, hi + pad)
    top, bot = ax.get_ylim()[1], ax.get_ylim()[0]

    for flavor, x, y in curves:
        st = STYLE[flavor]
        lbl = st["label"] if flavor not in seen else None
        seen.add(flavor)
        if flavor == "infinite":
            # draw the diverging path ONLY up to where it leaves the finite window (so its post-blowup
            # tail does not scribble a vertical line at t=0); mark the cutoff with an arrowhead pointing
            # the way it is going and the infinity symbol just inside the axis.
            outside = np.where((y > top) | (y < bot))[0]
            i = int(outside[0]) if len(outside) else len(x) - 1
            ax.plot(x[:i + 1], y[:i + 1], ls=st["ls"], color=st["color"], lw=1.8, alpha=0.9,
                    label=lbl, zorder=2)
            ax.scatter([x[0]], [y[0]], facecolor="gold", edgecolor="black", s=70, zorder=4, marker="o")
            if len(outside):
                # the line is simply cut off where it leaves the window (no end marker); the infinity
                # symbol just inside the axis says where it was heading.
                up = y[i] > top
                edge = top if up else bot
                dy = y[i] - y[i - 1] if i > 0 else 1.0       # crossing point, interpolated
                xe = x[i - 1] + (edge - y[i - 1]) / dy * (x[i] - x[i - 1]) if i > 0 and dy != 0 else x[i]
                inset = 0.05 * (top - bot)
                ax.text(xe, edge - inset if up else edge + inset, r"$\infty$", color=st["color"],
                        fontsize=13, ha="center", va="top" if up else "bottom", zorder=6)
        else:
            ax.plot(x, y, ls=st["ls"], color=st["color"], lw=1.8, alpha=0.9, label=lbl, zorder=2)
            ax.scatter([x[0]], [y[0]], facecolor="gold", edgecolor="black", s=70, zorder=4, marker="o")
            ax.scatter([x[-1]], [y[-1]], color=st["color"], edgecolor="black",
                       s=130, zorder=5, marker=st["marker"])               # endpoint (by flavor)

    ax.axvline(0.0, color="0.6", lw=1.5, zorder=1)          # the target t = 0
    ax.axvline(1.0, color="0.6", lw=1.5, zorder=1)          # the start system t = 1
    # a bullseye 🎯 just above t=0 marks the target system f(z)=0 (matplotlib's default fonts can't
    # render the emoji glyph, so we draw the dartboard the cartoon itself uses).
    tx = ax.get_xaxis_transform()                            # x in data coords, y in axes coords
    for size, col in [(380, "#c62828"), (190, "white"), (60, "#c62828")]:
        ax.scatter([0.0], [1.07], s=size, c=col, transform=tx, clip_on=False, zorder=10)
    ax.axvline(ENDGAME_BOUNDARY, color="purple", ls="--", lw=1.4)
    ax.text(ENDGAME_BOUNDARY, bot, " endgame boundary",
            color="purple", va="bottom", ha="left", fontsize=9, rotation=90)

    ax.set_xlim(-0.03, 1.03)          # t = 0 on the LEFT, t = 1 on the RIGHT (as in the cartoon)
    ax.set_xlabel(r"path variable $t$")
    ax.set_ylabel(r"generic real projection of the dehomogenized point")
    ax.set_title("The classic homotopy-continuation cartoon, from real data\n"
                 r"$f_1=(xy-3x+2)(x-4),\; f_2=y-x^2$")
    ax.scatter([], [], facecolor="gold", edgecolor="black", marker="o", s=70,
               label="start point ($t=1$)")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    ax.set_yticks([])                       # the height is a schematic projection; the scale is not meaningful
    ax.set_xticks([0.0, 1.0])               # only the two endpoints of the homotopy matter
    ax.set_xticklabels(["0", "1"])
    ax.grid(False)

    fig.tight_layout()
    svg, png = out_stem + ".svg", out_stem + ".png"
    fig.savefig(svg);  fig.savefig(png, dpi=150)
    plt.close(fig)
    return svg, png


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    paths = collect()
    from collections import Counter
    print("paths by endpoint flavor:", dict(Counter(f for f, _ in paths)))
    svg, png = make_plot(paths, os.path.join(out_dir, "classic_continuation_cartoon"))
    print("wrote", svg, "\n      ", png)


if __name__ == "__main__":
    main()
