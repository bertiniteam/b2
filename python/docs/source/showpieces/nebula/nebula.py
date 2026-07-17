"""The Nebula.

A showpiece: **one solve, every step**.  Where Homotopy Basins photographs *parameter* space --
many solves, one endpoint each -- the Nebula photographs *solution space in transit*.  Take a single
polynomial system with a lot of paths (cyclic-7: 5040 total-degree paths), track every one of them
from t=1 to t=0, keep **every accepted tracker step** along the way, project each C^n sample down to
a plane, and expose all of it onto one additive buffer.  Thousands of translucent trails pile up
into a gas cloud: bright knots where paths converge on solutions, filaments where bundles of paths
shear past the discriminant, voids where nothing travels.

It is a long exposure, and the exposure is literal:

    brightness  =  time spent

Each segment between two tracked samples carries the path-time ``|dt|`` the tracker spent crossing
it, spread evenly along the segment's footprint.  A path that crawls burns bright; a path that
races to infinity leaves a faint comet streak.  Brightness is *dwell*, spatially resolved -- and
because the energy of a segment is fixed regardless of how many samples the tracker chose to leave
on it, the picture is of the **path**, not of the stepper.  Nothing here is a synthetic texture:
every photon is real tracked data.

Channels (all honest)
---------------------
* **position**  -- a fixed complex-linear functional ``l(x) = sum a_k x_k``, plotted as
  ``(Re l, Im l)``.  Choosing ``l`` is this piece's composition knob (gamma was Basins').
* **brightness** -- dwell, as above.
* **colour**    -- temperature by ``log|t|``: cool where the path starts, white-hot in the endgame.
  Flow direction without arrows, and the knots at the solutions come out incandescent.
  (Note that colouring by ``arg l(x)`` would be *degenerate*: position already IS ``(Re l, Im l)``,
  so ``arg l`` is nothing but the screen polar angle -- a pinwheel carrying no information.)

Why the power-series endgame
----------------------------
The Cauchy endgame samples in *circles* around t=0, so ``|t|`` is deliberately multi-sheeted and
non-monotone, and a path's data can teleport back to the endgame boundary.  The power-series endgame
descends radially instead: measured over 5880 cyclic-5/6/7 paths, ``|t|`` is monotone on **every**
path and there are **zero** teleports.  The same lesson was learned the hard way by the Flight
Recorder showpiece (commit 0cbc7567).  It also keeps the piece cheap: under Cauchy, samples/path
explode superlinearly (259 -> 970 -> 11861 for n=5,6,7) and cyclic-7 becomes a ~29-hour render;
under power-series they stay flat (202 -> 462 -> 335) and cyclic-7 is **68 seconds**.

Determinism
-----------
There is a genuine random draw here (the total-degree start system and its gamma), unlike Basins
where every coefficient is exact.  So the seed IS the parameter: ``bertini.random.set_random_seed``
immediately before building, and ``num_threads = 1``.  Serial is not a compromise -- with a Python
observer attached, threading is ~7x *slower* (12 threads convoy on the GIL, one acquire per step).

Regenerated through ``tools/refresh_doc_artifacts.py`` (see that tool).  This is a raster showpiece:
PNG only (the image IS an accumulation buffer; there is no meaningful vector form), so it is exempt
from the tutorial figures' png+svg rule.  It is NOT a doctest -- the docs page embeds the
pre-rendered image.

Run standalone:  python nebula.py
"""

import argparse
import datetime
import os
import time

import numpy as np

import bertini
import bertini.nag_algorithm
from bertini import ZeroDimSolver, SolutionPathCollector

_OUT = os.path.dirname(os.path.abspath(__file__))

# supersampling: expose at SS x the output size, box-downsample for silky filaments
_SS = 2

# a segment whose |t| RATIO to the previous sample exceeds this is a teleport, not a step.
# Power-series never produces one (measured: 0 in 5880 paths), but the Cauchy endgame does -- a jump
# back to the endgame boundary is ~512x, while a legitimate Cauchy chord dips only ~1.95x.  |dt|
# alone CANNOT discriminate: a teleport's |dt| is ~0.0998, inside the max_step_size=0.1 cap.
_TELEPORT_RATIO = 3.0

# subsample spacing along a segment, in pixels, when rasterising the polyline
_STEP_PX = 0.5


# --- the subjects ---------------------------------------------------------------------------------

def cyclic_system(n):
    """The cyclic-n system: n! total-degree paths (cyclic-7 -> 5040).  Highly symmetric."""
    x = [bertini.Variable('x' + str(i)) for i in range(n)]
    y = list(x) + list(x)
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(x))
    for ii in range(n - 1):
        sys_.add_function(np.sum([np.prod(y[jj:jj + ii + 1]) for jj in range(n)]))
    sys_.add_function(np.prod(x) - 1)
    return sys_


def noon_system(n, c=None):
    """The noon-n benchmark -- Noonburg's neural network, ``3**n`` paths::

        x_i * (sum_{j != i} x_j^2)  -  c * x_i  +  1  =  0,     i = 1 .. n

    A Lotka-Volterra-style model of n interacting neurons (Noonburg, 1989), and a long-standing
    benchmark family for polynomial system solvers; the standard member is ``c = 11/10``, which is
    the default here.  Exact rational, per the coercion doctrine.
    """
    from fractions import Fraction
    c = Fraction(11, 10) if c is None else Fraction(c)
    x = [bertini.Variable('x' + str(i)) for i in range(n)]
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(x))
    cc = bertini.coefficient(c)
    for i in range(n):
        sys_.add_function(x[i] * np.sum([x[j]**2 for j in range(n) if j != i]) - cc * x[i] + 1)
    return sys_


def noonburg_like_system(n):
    """A Noonburg-*style* neural network, ``3**n`` paths -- and NOT the noon-n benchmark::

        x_i * (1 + sum_{j != i} x_j^2)  -  1  =  0,     i = 1 .. n

    Be precise about what this is, because it is the Nebula's teaching frame and it would be easy
    to mis-cite.  It is a sibling of :func:`noon_system`, not a member of it: writing noon's form as
    ``x_i sum_{j != i} x_j^2 - c x_i + 1``, this has the ``x_i`` coefficient of ``c = -1`` but the
    opposite constant.  Substituting ``x -> -x`` does carry it onto noon at ``c = -1``, so the two
    have mirror-image *solution sets*.

    Their **paths differ regardless**, which is what matters to a picture of paths: a total-degree
    homotopy runs from a FIXED start system to the target, so reflecting the target does not reflect
    the journey.  Measured, this system leaves 108k tracked samples where noon at ``c = -1`` leaves
    25k, and the two render quite differently.  This one fans into two broad wings; it was chosen by
    looking.
    """
    x = [bertini.Variable('x' + str(i)) for i in range(n)]
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(x))
    for i in range(n):
        sys_.add_function(x[i] * (1 + np.sum([x[j]**2 for j in range(n) if j != i])) - 1)
    return sys_


def katsura_system(n):
    """The Katsura-n system: 2^n paths.  A classic, and not symmetric."""
    x = [bertini.Variable('x' + str(i)) for i in range(n + 1)]

    def u(k):
        return x[abs(k)] if abs(k) <= n else None

    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(x))
    for i in range(n):
        terms = [u(j) * u(i - j) for j in range(-n, n + 1)
                 if u(j) is not None and u(i - j) is not None]
        sys_.add_function(np.sum(terms) - u(i))
    sys_.add_function(np.sum([u(j) for j in range(-n, n + 1) if u(j) is not None]) - 1)
    return sys_


def kuramoto_system(n, seed=0):
    """Equilibria of the Kuramoto model: ``n`` coupled oscillators, ``4**(n-1)`` paths.

    The Kuramoto model is how a crowd of oscillators with different natural frequencies pulls
    itself into synchrony -- fireflies flashing together, a power grid locking to 50 Hz.  Its
    equilibria satisfy

        omega_i  -  (K/n) * sum_j sin(theta_i - theta_j)  =  0 ,

    which is not polynomial until you set ``s_i = sin theta_i``, ``c_i = cos theta_i`` and carry
    ``s_i^2 + c_i^2 = 1`` along; then ``sin(theta_i - theta_j) = s_i c_j - c_i s_j`` and the whole
    thing is a square polynomial system.  That algebraization is exactly how the model's equilibria
    have been *counted* by homotopy continuation, so this is a system Bertini is genuinely for.

    Two constraints the mathematics imposes, both load-bearing:

    * **Gauge.** Rotating every phase by the same angle maps equilibria to equilibria, so the
      solution set is positive-dimensional (curves, not points) and a zero-dim solve would be
      meaningless.  Pin the last oscillator at ``theta = 0`` (``s = 0``, ``c = 1``) to quotient the
      symmetry out.
    * **Frequencies must sum to zero**, or no equilibrium exists at all: summing the equations
      kills the coupling term (it is antisymmetric in i, j) and leaves ``sum omega_i = 0``.  So the
      last frequency is *derived*, never drawn.

    The frequencies are exact rationals from a seeded RNG -- generic enough to break every symmetry
    of the system, and exact per the library's coercion doctrine.
    """
    from fractions import Fraction
    rng = np.random.default_rng(seed)
    omega = [Fraction(int(rng.integers(-9, 10)), 10) for _ in range(n - 1)]
    omega.append(-sum(omega))                    # sum omega_i = 0, or there is no equilibrium

    s = [bertini.Variable('s' + str(i)) for i in range(n - 1)]
    c = [bertini.Variable('c' + str(i)) for i in range(n - 1)]
    s_all = list(s) + [0]                        # the gauge: theta_{n-1} = 0
    c_all = list(c) + [1]

    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(s + c))
    coupling = bertini.coefficient(Fraction(1, n))          # K = 1
    for i in range(n - 1):
        interaction = None
        for j in range(n):
            term = s_all[i] * c_all[j] - c_all[i] * s_all[j]
            interaction = term if interaction is None else interaction + term
        sys_.add_function(bertini.coefficient(omega[i]) - coupling * interaction)
    for i in range(n - 1):
        sys_.add_function(s[i]**2 + c[i]**2 - 1)
    return sys_


def stewart_gough_system(seed=0):
    """Forward kinematics of a Stewart-Gough platform: 8192 paths, the classic 40 poses.

    This is the problem numerical algebraic geometry grew up on -- Bertini's own authors wrote the
    book on polynomial systems arising in kinematics.  A Stewart-Gough platform is a rigid plate
    held over a base by six legs of adjustable length (flight simulators, machine tools, telescope
    mounts).  *Inverse* kinematics is trivial: pick a pose, measure the legs.  **Forward**
    kinematics is the hard one -- given the six leg lengths, where is the platform?  The famous
    answer is that a general platform admits **40** poses for the same six lengths, and that count
    was settled by homotopy continuation.

    Unknowns: a unit quaternion ``q`` for the orientation and a translation ``t``, so seven of them.
    Equations: each leg pins a distance, ``|R(q) a_i + t - b_i|^2 = L_i^2``, plus ``|q|^2 = 1``.
    ``R(q)`` is quadratic in ``q``, so each leg equation has total degree 4 and the total-degree
    homotopy tracks ``4**6 * 2 = 8192`` paths for 40 answers -- a ratio that is itself the reason
    the field invented better start systems.

    The leg lengths are not invented: a pose is *chosen* and the lengths are computed from it
    exactly, so the system is guaranteed to have that pose among its solutions.  The chosen
    quaternion is an exact unit quaternion ``(1, 2, 2, 4)/5`` -- a Pythagorean quadruple, since an
    exact rational point on the unit sphere cannot be had by rounding.
    """
    from fractions import Fraction
    rng = np.random.default_rng(seed)

    def frac():
        return Fraction(int(rng.integers(-9, 10)), 5)

    a = [[frac() for _ in range(3)] for _ in range(6)]      # platform joints, in the moving frame
    b = [[frac() for _ in range(3)] for _ in range(6)]      # base joints, in the fixed frame

    # the pose we will hide in the answer: an EXACT unit quaternion (1,2,2,4)/5, |q|^2 = 25/25 = 1
    q_star = [Fraction(1, 5), Fraction(2, 5), Fraction(2, 5), Fraction(4, 5)]
    t_star = [frac() for _ in range(3)]

    def rot(q):
        """R(q) for a UNIT quaternion q, entries quadratic in q.  Works on Fractions or on nodes."""
        w, i, j, k = q
        return [[w*w + i*i - j*j - k*k, 2*(i*j - w*k), 2*(i*k + w*j)],
                [2*(i*j + w*k), w*w - i*i + j*j - k*k, 2*(j*k - w*i)],
                [2*(i*k - w*j), 2*(j*k + w*i), w*w - i*i - j*j + k*k]]

    # the leg lengths, computed EXACTLY from the chosen pose -- so it is genuinely a solution
    R_star = rot(q_star)
    leg_sq = []
    for ai, bi in zip(a, b):
        d = [sum(R_star[r][c] * ai[c] for c in range(3)) + t_star[r] - bi[r] for r in range(3)]
        leg_sq.append(sum(x * x for x in d))

    qv = [bertini.Variable(n) for n in ('q0', 'q1', 'q2', 'q3')]
    tv = [bertini.Variable(n) for n in ('t0', 't1', 't2')]
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(qv + tv))

    R = rot(qv)
    for ai, bi, Lsq in zip(a, b, leg_sq):
        f = None
        for r in range(3):
            d = R[r][0] * bertini.coefficient(ai[0]) \
                + R[r][1] * bertini.coefficient(ai[1]) \
                + R[r][2] * bertini.coefficient(ai[2]) \
                + tv[r] - bertini.coefficient(bi[r])
            f = d * d if f is None else f + d * d
        sys_.add_function(f - bertini.coefficient(Lsq))
    sys_.add_function(qv[0]**2 + qv[1]**2 + qv[2]**2 + qv[3]**2 - 1)
    return sys_


def mass_action_system(species, reactions, conserved, seed=0):
    """Steady states of a chemical reaction network under mass-action kinetics.

    The polynomials are *derived*, not transcribed: mass action is a fixed rule, so a network
    determines its own system.  For each reaction ``sum_i a_i X_i -> sum_i b_i X_i`` with rate
    constant ``k``, the flux is ``k * prod_i x_i**a_i``, and each species accumulates
    ``dx_i/dt = sum_reactions (b_i - a_i) * flux``.  Steady state sets every ``dx_i/dt = 0``.

    Those equations are never independent: a reaction network conserves things (total enzyme, total
    substrate), so the stoichiometric matrix is rank-deficient and the naive system is singular --
    a positive-dimensional solution set, useless to a zero-dim solver.  For each conservation law
    one ODE must be *replaced* by the corresponding linear conservation equation.  That is what
    ``conserved`` does, and it is what makes the system square and zero-dimensional.

    Parameters
    ----------
    species : list of str
        Species names, in order; they become the variables.
    reactions : list of (dict, dict)
        Each reaction as ``(reactants, products)``, mapping species name -> stoichiometry.  Rate
        constants are drawn as exact rationals from ``seed``.
    conserved : list of (list of str, str)
        Each conservation law as ``(species in the sum, the species whose ODE it replaces)``.  The
        conserved total is drawn as an exact rational.
    seed : int
        Seeds the rate constants and conserved totals.  Exact rationals throughout, per the
        library's coercion doctrine.
    """
    from fractions import Fraction
    rng = np.random.default_rng(seed)
    x = {s: bertini.Variable(s) for s in species}

    def rational():
        return bertini.coefficient(Fraction(int(rng.integers(1, 10)), 10))

    # flux of each reaction: k * prod reactant^stoichiometry
    fluxes = []
    for reactants, _ in reactions:
        flux = rational()
        for s, a in reactants.items():
            for _ in range(a):
                flux = flux * x[s]
        fluxes.append(flux)

    # dx_i/dt = sum_j (b_ij - a_ij) * flux_j
    ode = {s: None for s in species}
    for (reactants, products), flux in zip(reactions, fluxes):
        for s in set(reactants) | set(products):
            net = products.get(s, 0) - reactants.get(s, 0)
            if net == 0:
                continue
            term = net * flux
            ode[s] = term if ode[s] is None else ode[s] + term

    replaced = {victim for _, victim in conserved}
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup([x[s] for s in species]))
    for s in species:
        if s not in replaced and ode[s] is not None:
            sys_.add_function(ode[s])
    for members, _ in conserved:
        total = None
        for s in members:
            total = x[s] if total is None else total + x[s]
        sys_.add_function(total - rational())
    return sys_


def phosphorylation_network(n, seed=0):
    """The n-site distributive sequential phosphorylation network: ``2**(3n)`` paths.

    The workhorse of chemical reaction network theory, and the standard example of
    *multistationarity* -- a cell switch.  A kinase E walks a substrate up the ladder
    ``S_0 -> S_1 -> ... -> S_n`` one phosphate at a time, a phosphatase F walks it back down, and
    each step goes through an enzyme-substrate complex::

        S_i + E <-> ES_i -> S_{i+1} + E        (i = 0 .. n-1)
        S_i + F <-> FS_i -> S_{i-1} + F        (i = 1 .. n)

    Three things are conserved -- total kinase, total phosphatase, total substrate -- so three of
    the ODEs are redundant and are replaced by those conservation laws (see
    :func:`mass_action_system`).  n=3 gives 512 paths, n=4 gives 4096.
    """
    S = ['S{}'.format(i) for i in range(n + 1)]
    ES = ['ES{}'.format(i) for i in range(n)]
    FS = ['FS{}'.format(i + 1) for i in range(n)]
    species = S + ['E', 'F'] + ES + FS

    reactions = []
    for i in range(n):                                  # kinase: S_i + E <-> ES_i -> S_{i+1} + E
        reactions.append(({S[i]: 1, 'E': 1}, {ES[i]: 1}))
        reactions.append(({ES[i]: 1}, {S[i]: 1, 'E': 1}))
        reactions.append(({ES[i]: 1}, {S[i + 1]: 1, 'E': 1}))
    for i in range(1, n + 1):                           # phosphatase: S_i + F <-> FS_i -> S_{i-1} + F
        reactions.append(({S[i]: 1, 'F': 1}, {FS[i - 1]: 1}))
        reactions.append(({FS[i - 1]: 1}, {S[i]: 1, 'F': 1}))
        reactions.append(({FS[i - 1]: 1}, {S[i - 1]: 1, 'F': 1}))

    conserved = [(['E'] + ES, 'E'),                     # total kinase, free + bound
                 (['F'] + FS, 'F'),                     # total phosphatase
                 (S + ES + FS, S[0])]                   # total substrate, in every form
    return mass_action_system(species, reactions, conserved, seed=seed)


def dense_system(num_vars, degree, seed=0):
    """A random dense system: ``degree**num_vars`` paths, no structure and no symmetry at all.

    The note's own suggestion for the true nebula -- a system with nothing special about it.  Every
    coefficient is an EXACT rational drawn from a seeded RNG, per the library's coercion doctrine,
    so the system is reproducible and no float ever reaches the function tree.
    """
    from fractions import Fraction
    rng = np.random.default_rng(seed)
    x = [bertini.Variable('x' + str(i)) for i in range(num_vars)]
    sys_ = bertini.System()
    sys_.add_variable_group(bertini.VariableGroup(x))

    def monomials(d):
        """Every exponent tuple of total degree <= d."""
        out = [()]
        for _ in range(num_vars):
            out = [e + (k,) for e in out for k in range(d + 1) if sum(e) + k <= d]
        return out

    monos = monomials(degree)
    for _ in range(num_vars):
        f = None
        for e in monos:
            if sum(e) > degree:
                continue
            c = bertini.coefficient(Fraction(int(rng.integers(-9, 10)), 10)) \
                + bertini.I * bertini.coefficient(Fraction(int(rng.integers(-9, 10)), 10))
            term = c
            for xi, k in zip(x, e):
                if k:
                    term = term * xi**k
            f = term if f is None else f + term
        sys_.add_function(f)
    return sys_


_SYSTEMS = {
    'cyclic5': (lambda: cyclic_system(5), 'cyclic-5', 120),
    'cyclic6': (lambda: cyclic_system(6), 'cyclic-6', 720),
    'cyclic7': (lambda: cyclic_system(7), 'cyclic-7', 5040),
    'noon6': (lambda: noonburg_like_system(6), 'Noonburg-style, 6 neurons', 729),
    'noon6_benchmark': (lambda: noon_system(6), 'noon-6 benchmark (c = 11/10)', 729),
    'katsura8': (lambda: katsura_system(8), 'katsura-8', 256),
    'kuramoto5': (lambda: kuramoto_system(5), 'Kuramoto, 5 oscillators', 256),
    'kuramoto6': (lambda: kuramoto_system(6), 'Kuramoto, 6 oscillators', 1024),
    'kuramoto7': (lambda: kuramoto_system(7), 'Kuramoto, 7 oscillators', 4096),
    'phos3': (lambda: phosphorylation_network(3), '3-site phosphorylation', 512),
    'phos4': (lambda: phosphorylation_network(4), '4-site phosphorylation', 4096),
    'stewart': (lambda: stewart_gough_system(), 'Stewart-Gough platform', 8192),
    'dense2': (lambda: dense_system(2, 24), 'dense random 2-var deg-24', 576),
    'dense3': (lambda: dense_system(3, 9), 'dense random 3-var deg-9', 729),
}


# --- the engine: one solve, every step ------------------------------------------------------------

class _Camera(SolutionPathCollector):
    """Streaming meta-observer: compact each finished path to numpy, then drop it.

    ``SolutionPathCollector`` spins up a ``PathDataCollector`` per path on ``PathStarted`` and
    harvests it on ``PathComplete`` -- capturing EVERY path the solver starts, whatever its
    outcome.  Left alone it holds all of them to the end; cyclic-7's 1.7M samples as Python lists
    of boxed complexes would be gigabytes.  So subclass it (the attach/detach dance and its
    re-entrancy rules are subtle and already correct) and drain ``series`` as it fills: memory
    stays O(one path).

    The splat deliberately does NOT happen here.  Observers re-acquire the GIL on every event, so
    any real work in ``Observe`` serialises; compacting is cheap, rendering is not.
    """

    def __init__(self):
        super().__init__()
        self.points = []        # list of (n_k, n_vars+1) complex64, homogeneous
        self.times = []         # list of (n_k,) complex128

    def Observe(self, event):
        super().Observe(event)                 # the base does all the tracker attach/detach work
        while self.series:
            c = self.series.pop()
            self.points.append(c.points().astype(np.complex64))
            self.times.append(c.times())


def track(system_name, seed=2, max_step_size=None, mptype='double', endgame='powerseries'):
    """Solve once and keep every accepted step of every path.

    Parameters
    ----------
    system_name : str
        A key of :data:`_SYSTEMS`.
    seed : int
        The random seed.  This is the piece's only entropy: the total-degree start system and its
        gamma are drawn from it.  Set immediately before building -- a previous solve advances the
        RNG state and would change the picture.
    max_step_size : str or None
        Cap on the tracker's step size, as an EXACT rational -- pass a string ("0.005"), never a
        float.  Smaller caps sample the same path more finely (the mathematics is unchanged, and
        the dwell weighting is invariant to sampling density); ``None`` keeps the default 0.1.
    mptype : str
        ``'double'`` (fast, for prototyping) or ``'adaptive'``.
    endgame : str
        ``'powerseries'`` -- see the module docstring on why not Cauchy.

    Returns
    -------
    dict
        The cache: ``points`` (n_samples, n_vars+1) complex64 homogeneous coordinates, ``times``
        (n_samples,) complex128, ``path_id`` (n_samples,) int32, plus provenance.
    """
    build, label, _ = _SYSTEMS[system_name]
    bertini.random.set_random_seed(seed)
    bertini.recording(False)                   # millions of throwaway samples: do not archive them
    solver = ZeroDimSolver(build(), mptype=mptype, endgame=endgame)
    cfg = solver.get_config(bertini.nag_algorithm.ZeroDimConfig)
    cfg.num_threads = 1                        # serial: deterministic, and faster than threading here
    solver.set_config(cfg)
    if max_step_size is not None:
        solver.get_tracker().configure(stepping={'max_step_size': max_step_size})

    cam = _Camera()
    solver.add_observer(cam)
    t0 = time.perf_counter()
    solver.solve()
    dt = time.perf_counter() - t0

    path_id = np.concatenate([np.full(len(t), i, dtype=np.int32)
                              for i, t in enumerate(cam.times)])
    cache = dict(points=np.concatenate(cam.points), times=np.concatenate(cam.times),
                 path_id=path_id, system=system_name, label=label, seed=seed,
                 mptype=mptype, endgame=endgame,
                 max_step_size=str(max_step_size), n_paths=len(cam.times),
                 n_solutions=len(solver.finite_solutions()), seconds=dt)
    print(f'  tracked {cache["n_paths"]} paths, {len(path_id):,} samples, '
          f'{cache["n_solutions"]} finite solutions, {dt:.1f}s')
    return cache


def save_cache(cache, path):
    """Write a tracked cache to an .npz.  Tracking is the expensive part; rendering is not."""
    np.savez_compressed(path, **cache)
    print('  cached', path)


def load_cache(path):
    """Read a cache written by :func:`save_cache`."""
    z = np.load(path, allow_pickle=False)
    return {k: (z[k] if z[k].ndim else z[k].item()) for k in z.files}


# --- projection: C^n -> the plane -----------------------------------------------------------------
#
# These coefficients are applied POST-HOC, to tracked data, at render time.  They are part of the
# camera, not of the homotopy -- so unlike every other showpiece's constants they need not be exact
# values.  No node ever sees them.

def ell_dft(n, j):
    """The j-th cyclic character, ``a_k = exp(2 pi i j k / n)``.

    For cyclic-n this is an eigenvector of the cyclic shift, so the *solution set* maps to itself
    under rotation by 2 pi / n.  The paths are not symmetric though -- the start system and gamma
    break it -- so expect a near-mandala rather than an exact one.
    """
    return np.exp(2j * np.pi * j * np.arange(n) / n)


def ell_random(n, seed=0):
    """Generic unit coefficients: a projection that respects no symmetry of the system."""
    rng = np.random.default_rng(seed)
    return np.exp(2j * np.pi * rng.random(n))


def ell_order(nvars):
    """The Kuramoto **order parameter**, ``r e^{i psi} = (1/N) sum_j e^{i theta_j}``.

    The one projection in this file that is not a choice.  In the ``(s, c)`` coordinates
    ``e^{i theta_j} = c_j + i s_j``, so the order parameter -- the physical measure of how
    synchronised the oscillators are, ``r = 1`` locked and ``r = 0`` incoherent -- is *exactly a
    complex-linear functional of the variables*.  Projecting along it means the frame IS
    synchronisation space, and every trail is a path's journey toward (or away from) sync.

    The pinned oscillator contributes a constant ``1/N``, which merely translates the picture; the
    window's centre absorbs it.
    """
    n_osc = nvars // 2 + 1                       # nvars = 2 * (n_osc - 1) after the gauge
    a = np.zeros(nvars, dtype=complex)
    a[:n_osc - 1] = 1j / n_osc                   # the s_j
    a[n_osc - 1:] = 1.0 / n_osc                  # the c_j
    return a


_PROJECTIONS = {
    'order': ell_order,
    'dft1': lambda n: ell_dft(n, 1),
    'dft2': lambda n: ell_dft(n, 2),
    'dft3': lambda n: ell_dft(n, 3),
    'coord0': lambda n: np.eye(n, dtype=complex)[0],
    'random0': lambda n: ell_random(n, 0),
    'random1': lambda n: ell_random(n, 1),
}


def project(cache, ell):
    """Dehomogenize and apply ``l(x)``.  Returns the complex screen coordinate per sample.

    The solver homogenizes and patches, so the cached points carry a homogenizing coordinate in
    column 0 and are *bounded*; dividing it out here is what lets divergers fly off frame instead
    of overflowing the cache.  Paths that reach infinity drive that coordinate to zero, which
    yields inf/nan -- masked at splat time, not patched over.
    """
    pts = cache['points']
    with np.errstate(divide='ignore', invalid='ignore'):
        aff = pts[:, 1:] / pts[:, 0][:, None]
    return aff @ ell


# --- the exposure ---------------------------------------------------------------------------------

def _segment_energy(cache):
    """Per-segment path-time ``|dt|``, and a mask of which segments are real.

    A segment joins samples k and k+1 of the SAME path.  Rejected: cross-path joins, zero-energy
    steps, and teleports -- a solver reuses one tracker for the main track and the endgame
    sub-tracks, so a path's data is a concatenation and can jump.  Power-series never jumps
    (measured), but guard anyway; the guard is free and the failure is spectacular (a full-frame
    streak carrying the largest energy in the image).
    """
    t = cache['times']
    pid = cache['path_id']
    at = np.abs(t)
    energy = np.abs(np.diff(t))
    same = pid[1:] == pid[:-1]
    not_teleport = at[1:] <= _TELEPORT_RATIO * np.maximum(at[:-1], 1e-300)
    return energy, same & not_teleport & (energy > 0)


_RAMP = np.array([[1.00, 1.00, 0.95],                     # white-hot: arriving at the solution
                  [1.00, 0.62, 0.20],                     # orange
                  [0.85, 0.25, 0.75],                     # violet
                  [0.35, 0.35, 1.00],                     # blue
                  [0.13, 0.60, 1.00]])                    # cool: t = 1, the start system


# log10|t| stops for the ramp above.  These are FIXED, not fitted to the data, and that is
# deliberate -- it is the one tuning decision this piece really has, and it was made by looking.
#
# Measured on cyclic-7: 95% of the *energy* lies in log|t| in [-1.9, 0] (the endgame carries 24% of
# the samples but only 8.4% of the dwell).  The tempting move is to spread the five stops across
# that bulk so the whole ramp gets used.  Do not: it makes the picture WORSE.  Additive blending
# sums the hues that overlap in a pixel, so spreading the bulk across white/orange/violet/blue
# blends thousands of crossing trails into beige mud.  Coherence beats coverage.  Keep the bulk of
# the dwell inside ONE cool family (here [-1.5, 0] -> violet..cyan, so the main track is coloured
# but never muddy) and spend the warm end on the endgame alone, which is spatially concentrated at
# the solutions and so stays vivid instead of averaging out.
_TEMP_STOPS = np.array([-5.0, -3.0, -1.5, -0.7, 0.0])


def _temperature(logt, stops=_TEMP_STOPS):
    """Colour by ``log|t|``: cool where paths start, white-hot in the endgame.

    Coherent under additive blending -- thousands of *hues* would sum to grey, but a temperature
    ramp sums to a hotter temperature, which is what a long exposure of converging paths should do.
    Below the first stop (the deep endgame) everything clips to white-hot, which is where the knots
    at the solutions come from.
    """
    x = np.clip(logt, stops[0], stops[-1])
    return np.stack([np.interp(x, stops, _RAMP[:, k]) for k in range(3)], axis=-1)


def _phase_hue(cache, ok, ell2, value=1.0):
    """Hue from ``arg l2(x)`` for a SECOND, independent projection ``l2``.

    Not ``arg l(x)``: position already IS ``(Re l, Im l)``, so ``arg l`` is nothing but the screen
    polar angle -- a pinwheel that carries no information at all.  An independent ``l2`` is a
    genuinely different direction in C^n, invisible in the projection, and it varies smoothly along
    a path, so trails come out as coherent bands of colour rather than noise.
    """
    import matplotlib.colors as mcolors
    z2 = project(cache, ell2)[:-1][ok]
    hue = (np.angle(z2) / (2 * np.pi)) % 1.0
    hsv = np.stack([hue, np.full_like(hue, 0.85), np.full_like(hue, value)], axis=-1)
    return mcolors.hsv_to_rgb(hsv)


def _clip_to_frame(u0, v0, u1, v1, w, h):
    """Liang-Barsky: clip each segment to the frame, returning the surviving [s0, s1] fractions.

    Essential, not cosmetic: a diverger's segment can be a million pixels long, and subdividing it
    before clipping would allocate millions of subsamples for a line that crosses a few. Energy is
    uniform per unit length, so a clipped segment keeps energy * (s1 - s0) -- the photons that left
    the frame are simply lost, as in any camera.
    """
    du, dv = u1 - u0, v1 - v0
    s0 = np.zeros_like(u0)
    s1 = np.ones_like(u0)
    for p, q in ((-du, u0 - 0), (du, w - u0), (-dv, v0 - 0), (dv, h - v0)):
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(p != 0, q / p, np.inf)
        enter = (p < 0)
        s0 = np.where(enter, np.maximum(s0, r), s0)
        s1 = np.where(p > 0, np.minimum(s1, r), s1)
        # p == 0 and q < 0: the segment is parallel to and outside this edge
        s1 = np.where((p == 0) & (q < 0), -1.0, s1)
    return s0, s1


def _deposit(buf, uu, vv, weight, rgb):
    """Bilinear (anti-aliased) point splat of weighted photons into an (h, w, 3) buffer."""
    h, w = buf.shape[:2]
    x, y = uu - 0.5, vv - 0.5
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    fx, fy = x - x0, y - y0
    for dx in (0, 1):
        for dy in (0, 1):
            xi, yi = x0 + dx, y0 + dy
            a = weight * (fx if dx else 1 - fx) * (fy if dy else 1 - fy)
            ok = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h) & (a > 0) & np.isfinite(a)
            if not ok.any():
                continue
            flat = yi[ok] * w + xi[ok]
            for ch in range(3):
                buf[..., ch] += np.bincount(flat, weights=a[ok] * rgb[ok, ch],
                                            minlength=h * w).reshape(h, w)


def _splat_segments(buf, u0, v0, du, dv, energy, rgb):
    """Subdivide each segment along its footprint and deposit its energy.

    ``m = max(1, ceil(L / step))`` subsamples each carry ``E/m``, so the sum is exactly ``E`` for
    any length -- long segments spread to ``E/L`` per pixel (dwell per unit length = 1/speed),
    sub-pixel segments dump all of ``E`` into one pixel with no ``1/L -> inf`` singularity.
    """
    L = np.hypot(du, dv)
    m = np.maximum(1, np.ceil(L / _STEP_PX)).astype(np.int64)
    seg = np.repeat(np.arange(len(m)), m)                       # which segment each subsample is on
    off = np.arange(int(m.sum())) - np.repeat(np.cumsum(m) - m, m)
    f = (off + 0.5) / m[seg]
    _deposit(buf, u0[seg] + du[seg] * f, v0[seg] + dv[seg] * f,
             energy[seg] / m[seg], rgb[seg])


def expose(cache, ell, center, halfwidth, w, h, color='temperature', ell2=None):
    """The long exposure: splat every tracked segment into an (h, w, 3) float buffer.

    Each segment deposits its path-time energy ``E = |dt|``, spread evenly over
    ``m = max(1, ceil(L / step))`` subsamples along its rasterised footprint.  The sum over the
    subsamples is exactly ``E`` for any length, so:

    * long segments  -> ``E/L`` per pixel = dwell per unit length = 1/speed.  A literal long
      exposure of a moving particle.
    * sub-pixel segments (the endgame) -> ``m = 1``, all of ``E`` in one pixel.  No ``1/L -> inf``
      singularity: the pixel is the sensor, and motion below it is not resolvable.

    Subdividing is line rasterisation -- the geometry we already asserted by drawing a polyline
    through the samples, exactly as the Loom draws ribbons between its samples.  It invents no
    tracked data.  The energy is conserved regardless, which is what ``--selftest`` checks.
    """
    z = project(cache, ell)
    energy, ok = _segment_energy(cache)

    halfheight = halfwidth * h / w
    u = (z.real - (center[0] - halfwidth)) / (2 * halfwidth) * w
    v = (z.imag - (center[1] - halfheight)) / (2 * halfheight) * h

    u0, u1 = u[:-1], u[1:]
    v0, v1 = v[:-1], v[1:]
    ok = ok & np.isfinite(u0) & np.isfinite(u1) & np.isfinite(v0) & np.isfinite(v1)

    s0, s1 = _clip_to_frame(u0, v0, u1, v1, w, h)
    ok = ok & (s1 > s0)
    if not ok.any():
        return np.zeros((h, w, 3)), 0.0

    u0, u1, v0, v1 = u0[ok], u1[ok], v0[ok], v1[ok]
    s0, s1, E = s0[ok], s1[ok], energy[ok]
    du, dv = u1 - u0, v1 - v0
    # the surviving piece of each segment, and the share of its energy that piece carries
    cu0, cv0 = u0 + du * s0, v0 + dv * s0
    cdu, cdv = du * (s1 - s0), dv * (s1 - s0)
    E = E * (s1 - s0)

    logt = np.log10(np.maximum(np.abs(cache['times'][:-1][ok]), 1e-300))
    if color == 'temperature':
        rgb = _temperature(logt)
    elif color == 'phase2':
        rgb = _phase_hue(cache, ok, ell2)
    elif color == 'both':
        # hue from the independent direction, brightened toward white as the endgame closes in
        hot = np.clip((logt - _TEMP_STOPS[2]) / (_TEMP_STOPS[0] - _TEMP_STOPS[2]), 0, 1)
        rgb = _phase_hue(cache, ok, ell2) * (1 - hot)[:, None] + hot[:, None]
    else:
        raise ValueError(f'unknown colour mode {color!r}')

    buf = np.zeros((h, w, 3))
    _splat_segments(buf, cu0, cv0, cdu, cdv, E, rgb)
    # what the buffer SHOULD total if every photon landed: each segment writes E into 3 channels
    # scaled by its colour, so the expectation carries the colour weights too.
    return buf, float((E * rgb.sum(axis=1)).sum())


# --- rendering ------------------------------------------------------------------------------------

def _downsample(img, factor):
    """Box-average downsample of an (h, w, 3) image by an integer factor."""
    h, w = img.shape[:2]
    return img[:h - h % factor, :w - w % factor]\
        .reshape(h // factor, factor, w // factor, factor, -1).mean(axis=(1, 3))


def tone_map(buf, exposure=4.0, gamma=1.8, ss=_SS):
    """The HDR curve.  The accumulation IS the bloom; no bloom pass is needed.

    Downsample in LINEAR light first, so ``exposure`` means the same thing at any resolution, then
    ``1 - exp(-k*b)`` -- the classic film response, which is why the long exposure is literal.  Each
    channel saturates at its own rate, so the hottest knots bleach to white on their own.
    """
    img = _downsample(buf, ss)
    lit = img[img > 0]
    scale = np.percentile(lit, 99.5) if lit.size else 1.0
    img = 1 - np.exp(-exposure * img / max(scale, 1e-30))
    return np.clip(img ** (1 / gamma), 0, 1)


def render(cache, ell, center, halfwidth, w, h, out_png, exposure=4.0, gamma=1.8):
    """Expose, tone-map, write.  Returns the image."""
    import matplotlib.image as mimage
    buf, _ = expose(cache, ell, center, halfwidth, w * _SS, h * _SS)
    img = tone_map(buf, exposure=exposure, gamma=gamma)
    mimage.imsave(out_png, img, origin='lower')
    print('  wrote', out_png)
    return img


def suggest_window(cache, ell, quantile=99.0):
    """A dwell-weighted window for a projection: where do the photons actually land?

    Weight by energy, or the divergers -- which travel furthest and dwell least -- drag the frame
    out to infinity.  Prints a suggestion to paste in as a literal constant; auto-windowing at
    render time would make a committed PNG depend on run-to-run data.
    """
    z = project(cache, ell)
    energy, ok = _segment_energy(cache)
    ok = ok & np.isfinite(z[:-1])
    zz, ee = z[:-1][ok], energy[ok]
    order = np.argsort(np.abs(zz))
    cw = np.cumsum(ee[order]) / ee.sum()
    r = np.abs(zz[order][np.searchsorted(cw, quantile / 100.0)])
    cen = np.average(zz[order][cw <= quantile / 100.0].real), \
        np.average(zz[order][cw <= quantile / 100.0].imag)
    return (round(float(cen[0]), 3), round(float(cen[1]), 3)), round(float(r), 3)


# --- prototyping ----------------------------------------------------------------------------------

_PROTO = os.path.expanduser('~/nebula_prototype_images')


def _stamp(name):
    """A datetime-stamped path in the prototype image record."""
    os.makedirs(_PROTO, exist_ok=True)
    ts = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    return os.path.join(_PROTO, f'{name}_{ts}.png')


def scout(system_name='cyclic5', seed=2, w=480, h=270, max_step_size=None):
    """Contact sheet: ONE solve, every projection.

    The tracked data does not depend on ``l`` -- that is the whole point of caching the track -- so
    a single solve renders every candidate projection.  Each panel is captioned and the sheet is
    stamped into the prototype record.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print(f'scouting {system_name}:')
    cache = track(system_name, seed=seed, max_step_size=max_step_size)
    n = cache['points'].shape[1] - 1

    keys = list(_PROJECTIONS)
    cols = 3
    rows = (len(keys) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 2.7 * rows), facecolor='#05060a')
    for ax, key in zip(np.ravel(axes), keys):
        ell = _PROJECTIONS[key](n)
        center, halfwidth = suggest_window(cache, ell)
        buf, _ = expose(cache, ell, center, halfwidth, w * _SS, h * _SS)
        ax.imshow(tone_map(buf), origin='lower')
        ax.set_title(f'{key}   center={center} halfwidth={halfwidth}',
                     color='#c9d3e0', fontsize=8)
        print(f'  {key:8} center={center} halfwidth={halfwidth}')
    for ax in np.ravel(axes):
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f'{cache["label"]}  --  {cache["n_paths"]} paths, '
                 f'{len(cache["path_id"]):,} samples, seed {seed}', color='#e8eef7')
    out = _stamp(f'scout_{system_name}')
    fig.savefig(out, dpi=110, facecolor='#05060a', bbox_inches='tight')
    plt.close(fig)
    print('  wrote', out)
    return cache


def selftest(system_name='cyclic5'):
    """Assert the exposure conserves energy -- the sharpest invariant this piece has.

    Every photon in the buffer is some segment's path-time.  Two checks:

    1. **Synthetic, strictly interior**: conservation must be EXACT (to float tolerance).  This
       isolates the subdivision and the bilinear deposit from any real-data effect.  Segments of
       wildly different lengths -- including sub-pixel ones -- must all conserve.
    2. **Real tracked data**: the buffer will fall slightly short, because a bilinear splat at the
       frame border drops the neighbours that lie outside it.  That is honest photon loss (a camera
       does the same), so it is reported, not asserted away.
    """
    rng = np.random.default_rng(0)
    h, w = 64, 96
    n = 400
    u0 = rng.uniform(8, w - 8, n)
    v0 = rng.uniform(8, h - 8, n)
    ang = rng.uniform(0, 2 * np.pi, n)
    length = 10 ** rng.uniform(-1.5, 1.2, n)          # 0.03 px (sub-pixel) .. 16 px
    du, dv = length * np.cos(ang), length * np.sin(ang)
    E = rng.uniform(0.1, 5.0, n)
    rgb = rng.uniform(0.2, 1.0, (n, 3))
    buf = np.zeros((h, w, 3))
    _splat_segments(buf, u0, v0, du, dv, E, rgb)
    want = float((E * rgb.sum(axis=1)).sum())
    got = float(buf.sum())
    print(f'  synthetic: want {want:.9f}, got {got:.9f}, rel err {abs(got - want) / want:.2e}')
    assert abs(got - want) / want < 1e-9, 'the splat does not conserve energy'
    assert (m := (length < 1).sum()) > 0, 'test needs sub-pixel segments to be meaningful'
    print(f'             ({m} of {n} segments were sub-pixel)')

    cache = track(system_name, seed=2)
    ell = _PROJECTIONS['dft1'](cache['points'].shape[1] - 1)
    center, halfwidth = suggest_window(cache, ell)
    buf, expected = expose(cache, ell, center, halfwidth, 480, 270)
    got = float(buf.sum())
    print(f'  tracked  : want {expected:.6f}, got {got:.6f}, '
          f'{100 * (1 - got / expected):.3f}% lost at the frame border')
    assert got <= expected * (1 + 1e-9), 'the exposure created energy from nowhere'
    assert got > expected * 0.9, 'the exposure lost more than 10% -- that is not border loss'
    print('  PASS')


# --- the two frames -------------------------------------------------------------------------------
#
# The windows are LITERAL CONSTANTS, deliberately.  suggest_window() proposes them from the tracked
# data and prints them; they are then pasted here.  Auto-windowing at render time would make a
# committed PNG depend on run-to-run data and churn the diff on every regeneration.  Constants are
# reviewable.  Both frames project onto coord0 -- the first coordinate's own complex plane -- which
# was chosen by looking at contact sheets, not by argument.

def teaching_frame(out):
    """A Noonburg-style network on 6 neurons: 729 paths, 717 solutions, ~3s.  The lesson, legible.

    Few enough paths that individual strands stay separate, so you can see what the piece IS: each
    filament is one tracked path, fanning out of the start system and converging on a solution,
    brightening as it slows.  The wing-like plumes are bundles of paths that travel together.
    """
    print('teaching frame (noonburg-6):')
    cache = track('noon6', seed=2)
    ell = _PROJECTIONS['coord0'](cache['points'].shape[1] - 1)
    render(cache, ell, center=(0.1, -0.502), halfwidth=1.357, w=960, h=540,
           out_png=out, exposure=2.0, gamma=2.4)


def showpiece_frame(out):
    """Kuramoto, 7 oscillators: 4096 paths, 124 equilibria, ~40s.  The show-off.

    Dense enough to be a cloud rather than a set of strands.  Two blazing knots where whole bundles
    of paths converge, and thousands of filaments streaming between them: a solve, photographed.
    """
    print('showpiece frame (Kuramoto-7):')
    cache = track('kuramoto7', seed=2)
    ell = _PROJECTIONS['coord0'](cache['points'].shape[1] - 1)
    render(cache, ell, center=(0.278, 0.322), halfwidth=5.418, w=960, h=540,
           out_png=out, exposure=2.0, gamma=2.4)


def main():
    """Generate both frames next to this script, or scout/selftest while prototyping."""
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--scout', metavar='SYSTEM', nargs='?', const='cyclic5',
                   help='contact sheet of every projection, from ONE solve (prototyping)')
    p.add_argument('--selftest', action='store_true', help='assert the exposure conserves energy')
    p.add_argument('--seed', type=int, default=2)
    p.add_argument('--max-step-size', default=None,
                   help='exact rational as a string, e.g. "0.005"')
    args = p.parse_args()

    if args.selftest:
        selftest()
    elif args.scout:
        scout(args.scout, seed=args.seed, max_step_size=args.max_step_size)
    else:
        teaching_frame(os.path.join(_OUT, 'nebula_teaching.png'))
        showpiece_frame(os.path.join(_OUT, 'nebula.png'))


if __name__ == '__main__':
    main()
