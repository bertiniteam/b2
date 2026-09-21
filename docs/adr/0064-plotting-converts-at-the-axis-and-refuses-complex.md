# ADR-0064: plotting converts at the axis, and refuses complex data there

**Status:** Accepted
**Date:** 2026-09-21

## Context

A library whose numbers cannot be plotted does not get used.  Most of pyplot already worked on
`real_mp`, because the type converts to a float and matplotlib asks each element for one.  The
calls that compute against plain floats before they draw did not: `hist`, `bar` with a float
width, `boxplot` all failed with a numpy dtype promotion error, because `real_mp` and `float64`
have no common dtype.

That missing promotion is deliberate and stays missing.  Registering a safe cast to `float64`
would let numpy pick it as the common type in ordinary arithmetic, so `mp_array * 0.5` would
quietly become a double array and every digit past the sixteenth would be gone -- silently, in
computation, which is the opposite of what a multiprecision library is for.

Complex data was worse than unplottable.  numpy's `ndarray.real` / `.imag` are hardwired to its
own three built-in complex types (`PyArray_ISCOMPLEX`); on any other complex dtype `.real`
returns the array itself and `.imag` returns an array of zeros, with no error, no warning, and
no hook for a user dtype to correct or even detect it.  ADR-0006's companion guard
(`bertini._numpy_guard`) already makes `np.real` / `np.imag` / `np.angle` raise, because those
are plain Python functions we can wrap -- but the attribute spelling is C-level and untouchable.
So `scatter(points.real, points.imag)`, which is how everyone plots a complex array, drew every
point on the x axis and looked like an answer.  Measured, before this change:

    offsets drawn: [[0.142, 0.0], [0.285, 0.0], [0.428, 0.0], ...]
    truth:         [(0.142, 1.0), (0.285, 2.0), (0.428, 3.0), ...]

## Decision

**Conversion happens at the axis, through matplotlib's units registry, and nowhere else.**
`bertini._matplotlib_bridge` registers a `ConversionInterface` for `real_mp`, which hands
matplotlib `float64` for data bound for an axis.  A converter is consulted only for axis data,
which is the one place where dropping digits is right: a screen resolves about three of them,
and the values themselves are never touched.  Image data, contour levels and marker sizes do not
pass that way and are cast explicitly with `.astype(float)`.

**Complex data is refused at the axis rather than converted.**  The converter registered for
`complex_mp` raises a `TypeError` naming `bertini.real` / `bertini.imag`, and the real converter
refuses complex input too, because an axis keeps the first converter it is given.  This works as
a guard precisely because numpy's `.real` of a `complex_mp` array is still `complex_mp`: the
wrong spelling arrives at the axis still carrying the complex dtype, so refusing complex catches
exactly the mistake that would otherwise be silent.  A `Solution` is unaffected -- it overrides
`.real` / `.imag` correctly, and hands the axis real data.

**Registration is automatic and does not import matplotlib.**  `install()` registers immediately
if `matplotlib.units` is already loaded, and otherwise leaves a `sys.meta_path` finder that wraps
that module's loader and registers as it finishes executing, then removes itself.  Import order
does not matter, and bertini still imports on a machine with no matplotlib at all.

## Consequences

- **Do not "fix" the complex converter to plot real parts.**  Drawing the real parts of complex
  data is the failure mode this exists to prevent, not a convenience to restore.  If a caller
  wants them, `bertini.real(z)` says so out loud.
- **Do not register a safe or same-kind cast between the mp dtypes and `float64`** to make
  `imshow`, `contour`, `boxplot` or `violinplot` work.  Those refuse because
  `np.can_cast(real_mp, float, "same_kind")` is `False`, and it must stay `False`: the same
  registration that satisfies matplotlib would silently downcast ordinary arithmetic.
  `.astype(float)` at the call site is the supported spelling.
- Converting at the axis means tick labels and limits are float64.  That is display, and the
  stored values keep every digit.
- pandas screens columns by dtype before matplotlib is involved, so `DataFrame.plot()` and
  `.plot.hist()` still refuse an object column of mp values; `.plot.scatter` goes through the
  converter and works.  Libraries that serialize rather than draw (plotly) refuse mp values with
  their own error, which is loud and therefore fine.
