"""Homotopy Basins.

A showpiece: **every pixel is a zero-dimensional solve**.  Take the one-parameter family

    f(x; c) = x^d - d*x - c

and rasterize a window of the complex c-plane.  For each pixel, run the genuine total-degree
gamma-trick homotopy

    H(x, t) = gamma * t * (x^d - 1)  +  (1 - t) * (x^d - d*x - c)

and track all d start roots (the roots of unity) from t=1 to t=0 with the
``DoublePrecisionTracker``.  Nothing in the image is a synthetic texture -- every channel is
tracked data:

    * hue        -- the phase of the *start-end correlation* ``s(c) = sum_k zeta^k x_k``:
                    each path contributes (its start root) * (its landing point), with
                    zeta = exp(2 pi i / d).  The unweighted sum of landings would be
                    permutation-blind (here identically zero -- it is the x^(d-1) coefficient);
                    the start-root weights make s holomorphic in c inside each basin and make it
                    jump exactly where the homotopy's swept discriminant permutes which root
                    each path reaches;
    * brightness -- the tracker's own total step count.  The blazing arcs are the discriminant
                    of H swept through system space: the set of c for which the straight-line
                    homotopy passes through a singular system at some t.  The adaptive stepper
                    piles up tiny steps exactly there;
    * stripes    -- level bands of log|s| (honest domain coloring; cosmetic, not singularity);
                    the rainbow whirlpools are the zeros of s, coiled at the roots of the arcs;
    * white speckle -- pixels where a path actually failed near a singular system.

The geometry is the gamma trick, photographed.  The family's branch points form the ring
c = -(d-1) * zeta with zeta^(d-1) = 1; each grows one glowing arc, and every arc trails off to
infinity in the -gamma direction (rotating gamma literally re-aims the comet tails).

Four frames are produced (each window contains the family's whole branch-point ring; each
degree is rendered twice from one computation, full-bleed and as an annotated map with labeled
Re(c)/Im(c) axes and every arc center c = -(d-1) zeta marked on its ring):

    * homotopy_basins_teaching.png / homotopy_basins_teaching_annotated.png
                                    -- d = 5: four branch points, four arcs, the lesson legible.
    * homotopy_basins.png / homotopy_basins_annotated.png
                                    -- d = 9: the eight-arc comet cluster, the show-off.

Every coefficient fed to the function tree is EXACT (fractions.Fraction; pixel coordinates are
snapped to rationals), per the library's coercion doctrine.  There is no randomness anywhere --
gamma is a fixed exact unit -- so the tracked data is fully deterministic.

Regenerated through ``tools/refresh_doc_artifacts.py`` (see that tool).  This is a raster
showpiece: PNG only (the image IS a raster of solves; there is no meaningful vector form), so it
is exempt from the tutorial figures' png+svg rule.  It is NOT a doctest -- the docs page embeds
the pre-rendered image.  It is also the heaviest doc artifact: the two frames together are about
17 million tracked paths (roughly 10-15 minutes on 12 cores).

Run standalone:  python homotopy_basins.py
"""

import os
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np

_OUT = os.path.dirname(os.path.abspath(__file__))

# pixel coordinates are snapped onto this exact grid (coefficients must be exact values)
_DENOM = 10**6

# supersampling: compute at SS x the output size, box-downsample for silky bands and arcs
_SS = 2


# --- the engine: one total-degree homotopy per pixel --------------------------------------------

_worker = {}

def _worker_init(degree, gamma_re, gamma_im):
    """Per-process setup: import bertini once and remember the frame constants."""
    import bertini
    import bertini.tracking  # noqa: F401  (registers the tracking submodule)
    bertini.recording(False)              # millions of throwaway solves: do not archive them
    _worker['bertini'] = bertini
    _worker['degree'] = degree
    _worker['gamma'] = (gamma_re, gamma_im)
    _worker['starts'] = [np.exp(2j * np.pi * k / degree) for k in range(degree)]
    _worker['zeta'] = np.exp(2j * np.pi * np.arange(degree) / degree)


def _cnode(re_exact, im_exact):
    """An exact complex coefficient node re + i*im (both fractions.Fraction or int)."""
    bertini = _worker['bertini']
    return bertini.coefficient(re_exact) + bertini.I * bertini.coefficient(im_exact)


def _track_pixel(c_re, c_im):
    """All d total-degree paths for the target x^d - d*x - c at one exact pixel value c.

    Returns (start-end correlation s, total_steps, num_failures)."""
    bertini = _worker['bertini']
    tracking = bertini.tracking
    d = _worker['degree']

    x = bertini.Variable('x')
    t = bertini.Variable('t')
    H = bertini.System()
    H.add_function(_cnode(*_worker['gamma']) * t * (x**d - 1)
                   + (1 - t) * (x**d - d * x - _cnode(c_re, c_im)))
    H.add_path_variable(t)
    H.add_variable_group(bertini.VariableGroup([x]))

    tracker = bertini.DoublePrecisionTracker(H)
    tracker.setup(tracking.Predictor.RK4, 1e-6, 1e5,
                  tracking.SteppingConfig(), tracking.NewtonConfig())

    lands = np.zeros(d, dtype=complex)
    steps = 0
    fails = 0
    end = np.zeros(1, dtype=complex)
    for k, root in enumerate(_worker['starts']):
        code = tracker.track_path(end, complex(1, 0), complex(0, 0), np.array([root]))
        steps += tracker.num_total_steps_taken()
        if str(code) != 'Success':
            fails += 1
        lands[k] = end[0]
    return (lands * _worker['zeta']).sum(), steps, fails


def _row(args):
    """One raster row: track every pixel, reduce to (correlation s, steps, failures) arrays."""
    j, c_res, c_im = args
    finger = np.zeros(len(c_res), dtype=complex)
    steps = np.zeros(len(c_res), dtype=np.int32)
    fails = np.zeros(len(c_res), dtype=np.int8)
    for i, c_re in enumerate(c_res):
        finger[i], steps[i], fails[i] = _track_pixel(c_re, c_im)
    return j, finger, steps, fails


def compute(degree, gamma, center, halfwidth, w, h):
    """Track the whole window (w x h pixels, all d paths each).  Returns per-pixel arrays
    (correlation s, steps, failures)."""
    halfheight = halfwidth * h / w
    xs = [Fraction(round((center[0] - halfwidth + 2 * halfwidth * i / (w - 1)) * _DENOM), _DENOM)
          for i in range(w)]
    ys = [Fraction(round((center[1] - halfheight + 2 * halfheight * j / (h - 1)) * _DENOM), _DENOM)
          for j in range(h)]
    finger = np.zeros((h, w), dtype=complex)
    steps = np.zeros((h, w), dtype=np.int32)
    fails = np.zeros((h, w), dtype=np.int8)
    t0 = time.perf_counter()
    with Pool(os.cpu_count(), initializer=_worker_init,
              initargs=(degree, gamma[0], gamma[1])) as pool:
        jobs = [(j, xs, ys[j]) for j in range(h)]
        for n, (j, fg, st, fl) in enumerate(pool.imap_unordered(_row, jobs, chunksize=1)):
            finger[j] = fg; steps[j] = st; fails[j] = fl
            if n % max(1, h // 12) == 0:
                print(f'  row {n + 1}/{h}  ({time.perf_counter() - t0:.0f}s)', flush=True)
    print(f'  tracked {w}x{h} pixels x {degree} paths in {time.perf_counter() - t0:.0f}s')
    return finger, steps, fails


# --- rendering -----------------------------------------------------------------------------------

def _blur(channel, sigma):
    """FFT gaussian blur of one 2-D channel."""
    fy = np.fft.fftfreq(channel.shape[0])[:, None]
    fx = np.fft.fftfreq(channel.shape[1])[None, :]
    kernel = np.exp(-2 * (np.pi * sigma) ** 2 * (fx ** 2 + fy ** 2))
    return np.real(np.fft.ifft2(np.fft.fft2(channel) * kernel))


def _downsample(img, factor):
    """Box-average downsample of an (h, w, 3) image by an integer factor."""
    h, w = img.shape[:2]
    return img[:h - h % factor, :w - w % factor]\
        .reshape(h // factor, factor, w // factor, factor, -1).mean(axis=(1, 3))


def _shade(finger, steps, fails):
    """Map the tracked channels to an RGB image: hue from the phase of the start-end
    correlation s, brightness from tracker effort, stripes from |s|, speckle from failures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.colors as mcolors

    w = steps.shape[1]
    hue = (np.angle(finger) / (2 * np.pi)) % 1.0

    effort = np.log1p(steps.astype(float))
    effort = _blur(effort, sigma=max(0.7, w / 900))      # tame per-pixel step-count speckle
    lo, hi = np.percentile(effort, 5), np.percentile(effort, 99.7)
    stress = np.clip((effort - lo) / (hi - lo + 1e-9), 0, 1)

    mag = np.abs(finger)
    band = 0.5 + 0.5 * np.cos(2 * np.pi * 2.2 * np.log(mag + 1e-12))
    hue = (hue + 0.07 * np.log(mag + 1e-12)) % 1.0       # iridescent drift across the bands

    val = np.clip(0.08 + 0.22 * band + 0.30 * stress + 0.55 * stress ** 3, 0, 1)
    sat = np.clip(1.0 - 0.75 * stress ** 4, 0.35, 1.0)   # the hottest arcs bleach to white
    rgb = mcolors.hsv_to_rgb(np.stack([hue, sat, val], axis=-1))

    # bloom: the stress field re-added as soft glow in the local hue (screen blend)
    core = rgb * (stress ** 3)[..., None]
    glow = np.stack([_blur(core[..., k], sigma=max(1.5, w / 300)) for k in range(3)], axis=-1)
    img = 1 - (1 - rgb) * (1 - np.clip(1.3 * glow, 0, 1))
    img = np.clip(img + (fails > 0)[..., None] * 0.55, 0, 1)     # failure speckle, white-hot

    return np.clip(_downsample(img, _SS), 0, 1)


def render(finger, steps, fails, out_png):
    """The full-bleed frame: just the shaded image, no chrome."""
    import matplotlib.image as mimage
    mimage.imsave(out_png, _shade(finger, steps, fails), origin='lower')
    print('  wrote', out_png)


def render_annotated(finger, steps, fails, degree, extent, out_png):
    """The same tracked data as a *map*: labeled Re(c)/Im(c) axes, and every arc center --
    the family's branch points c = -(d-1) zeta, zeta^(d-1) = 1 -- marked on the ring."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    img = _shade(finger, steps, fails)
    bg = '#05060a'
    plt.rcParams.update({'figure.facecolor': bg, 'axes.facecolor': bg})
    fig, ax = plt.subplots(figsize=(12, 7.4))
    ax.imshow(img, origin='lower', extent=extent, interpolation='bilinear')

    radius = degree - 1
    branch = -radius * np.exp(2j * np.pi * np.arange(degree - 1) / (degree - 1))
    ax.scatter(branch.real, branch.imag, s=130, facecolors='none',
               edgecolors='white', linewidths=1.6, zorder=5)
    for b in branch:
        if abs(b.imag) < 1e-9:                              # on the real axis: -8 or 8
            label = f'${round(b.real)}$'
        elif abs(b.real) < 1e-9:                            # on the imaginary axis: -8i or 8i
            label = f'${round(b.imag)}i$'.replace('1i', 'i')
        else:                                               # diagonal: exact polar form
            label = f'${radius}e^{{{np.angle(b) / np.pi:+.2g}\\pi i}}$'
        ax.annotate(label, (b.real, b.imag), textcoords='offset points', xytext=(9, 9),
                    color='white', fontsize=9)
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot((degree - 1) * np.cos(theta), (degree - 1) * np.sin(theta),
            color='white', lw=0.6, ls=':', alpha=0.45,
            label=f'branch-point ring  $|c| = {degree - 1}$')

    ax.set_xlabel(r'$\operatorname{Re}\, c$', color='#c9d3e0')
    ax.set_ylabel(r'$\operatorname{Im}\, c$', color='#c9d3e0')
    ax.tick_params(colors='#8896aa')
    for spine in ax.spines.values():
        spine.set_color('#223')
    ax.set_title(f'arc centers: the branch points  $c = -{degree - 1}\\zeta$,  '
                 f'$\\zeta^{{{degree - 1}}} = 1$   (each arc trails toward $-\\gamma\\infty$)',
                 color='#e8eef7', fontsize=11)
    leg = ax.legend(loc='lower right', framealpha=0.15, labelcolor='#c9d3e0', fontsize=9)
    leg.get_frame().set_edgecolor('#223')
    fig.savefig(out_png, dpi=110, facecolor=bg, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', out_png)


# --- the two frames ------------------------------------------------------------------------------

def _extent(center, halfwidth, w, h):
    """The imshow extent [x0, x1, y0, y1] matching compute()'s pixel grid."""
    halfheight = halfwidth * h / w
    return [center[0] - halfwidth, center[0] + halfwidth,
            center[1] - halfheight, center[1] + halfheight]


def teaching_frame(out, out_annotated):
    """d = 5: four branch points (c = -4 zeta, zeta^4 = 1), four arcs -- the lesson legible.
    The window contains the whole branch-point ring |c| = 4; the same tracked data is rendered
    full-bleed and as an annotated map with the four arc centers marked."""
    print('teaching frame (d=5):')
    center, halfwidth = (1.2, -0.2), 9.0
    w, h = 800 * _SS, 450 * _SS
    finger, steps, fails = compute(degree=5,
                                   gamma=(Fraction(-24, 25), Fraction(7, 25)),
                                   center=center, halfwidth=halfwidth, w=w, h=h)
    render(finger, steps, fails, out)
    render_annotated(finger, steps, fails, 5, _extent(center, halfwidth, w, h), out_annotated)


def showpiece_frame(out, out_annotated):
    """d = 9: the eight-arc comet cluster streaming in the -gamma direction -- the show-off.
    The window contains the whole branch-point ring |c| = 8; the same tracked data is rendered
    twice, full-bleed and as an annotated map with the eight arc centers marked."""
    print('showpiece frame (d=9):')
    center, halfwidth = (2.5, -0.3), 16.4
    w, h = 1200 * _SS, 675 * _SS
    finger, steps, fails = compute(degree=9,
                                   gamma=(Fraction(-24, 25), Fraction(7, 25)),
                                   center=center, halfwidth=halfwidth, w=w, h=h)
    render(finger, steps, fails, out)
    render_annotated(finger, steps, fails, 9, _extent(center, halfwidth, w, h), out_annotated)


def main():
    """Generate all four frames next to this script."""
    teaching_frame(os.path.join(_OUT, 'homotopy_basins_teaching.png'),
                   os.path.join(_OUT, 'homotopy_basins_teaching_annotated.png'))
    showpiece_frame(os.path.join(_OUT, 'homotopy_basins.png'),
                    os.path.join(_OUT, 'homotopy_basins_annotated.png'))


if __name__ == '__main__':
    main()
