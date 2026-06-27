"""Test-bed polynomial systems for the external-solver comparison benchmark.

Every system is built in Python with pybertini and is genuinely zero-dimensional, so the
solution count is a known constant we can use as a correctness check (see ``expected_count``).
The driver emits each system *once* to a classic Bertini-1 input file via
``System.to_classic_input(...)`` and feeds that single file to every solver, so the problem and
all tracking settings are identical across solvers by construction.

Add a new family by writing a generator that returns a built ``bertini.System`` and registering
it in ``SYSTEMS``.  Keep generators free of any solver- or timing-specific concerns.
"""

import functools
import operator

import bertini as pb


def _product(factors):
    """Symbolic product of a non-empty list of nodes."""
    return functools.reduce(operator.mul, factors)


def _sum(terms):
    """Symbolic sum of a non-empty list of nodes."""
    return functools.reduce(operator.add, terms)


# --------------------------------------------------------------------------------------------
# Families
# --------------------------------------------------------------------------------------------

def cyclic(n):
    """The cyclic-n roots system in n variables.

    Equations: for k = 1 .. n-1 the sum of all length-k cyclic products of the variables, plus
    the closing relation (product of all variables) - 1.

    cyclic-n is zero-dimensional only when n is squarefree (Backelin); n in {5, 6, 7, 10, 11}
    are the small squarefree cases.  n in {4, 8, 9, 12, ...} have positive-dimensional
    components and are NOT valid zero-dim test cases -- don't register them.
    """
    x = pb.variables('x', n)               # x0 .. x_{n-1}
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(x))

    for k in range(1, n):
        terms = [_product([x[(i + j) % n] for j in range(k)]) for i in range(n)]
        sys.add_function(_sum(terms))
    sys.add_function(_product(x) - 1)
    return sys


def katsura(n):
    """The Katsura-n system in n+1 variables (x0 .. xn), with 2**n solutions.

    Uses the symmetric-index convention: x_{-i} = x_i, and x_i = 0 for |i| > n.  For
    m = 0 .. n-1 the convolution  sum_j x_{|j|} x_{|m-j|} - x_m = 0, plus the normalization
    x0 + 2*(x1 + ... + xn) - 1 = 0.
    """
    x = pb.variables('x', n + 1)           # x0 .. xn
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(x))

    for m in range(n):
        terms = []
        for j in range(-n, n + 1):
            a, b = abs(j), abs(m - j)
            if a <= n and b <= n:
                terms.append(x[a] * x[b])
        sys.add_function(_sum(terms) - x[m])
    sys.add_function(x[0] + 2 * _sum(x[1:]) - 1)
    return sys


def diagonal(n, degree=3):
    """The diagonal warmup family x_i**degree - c_i, with degree**n solutions.

    Trivially zero-dimensional with analytically known roots -- a fast sanity tier mirroring the
    existing benchmark/inputs/*.b2 systems.  The c_i are distinct small primes.
    """
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n > len(primes):
        raise ValueError(f"diagonal: need {n} distinct constants, have {len(primes)} primes")
    x = pb.variables('x', n)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(x))
    for i in range(n):
        sys.add_function(x[i] ** degree - primes[i])
    return sys


# --------------------------------------------------------------------------------------------
# Registry: name -> (builder, expected solution count)
# --------------------------------------------------------------------------------------------
# expected_count is the number of isolated complex solutions of the (zero-dimensional) system;
# it is the ground truth for the correctness / agreement check, never part of any timing.

SYSTEMS = {
    # diagonal warmup tier (fast sanity)
    "diag3":    (lambda: diagonal(3), 27),
    "diag5":    (lambda: diagonal(5), 243),
    "diag6":    (lambda: diagonal(6), 729),

    # cyclic-n (squarefree n only: genuinely zero-dimensional)
    "cyclic5":  (lambda: cyclic(5), 70),
    "cyclic6":  (lambda: cyclic(6), 156),
    "cyclic7":  (lambda: cyclic(7), 924),

    # Katsura-n (2**n solutions)
    "katsura3": (lambda: katsura(3), 8),
    "katsura4": (lambda: katsura(4), 16),
    "katsura5": (lambda: katsura(5), 32),
    "katsura6": (lambda: katsura(6), 64),
}

# Sensible default subset for a quick run: one of each family, mid-sized.
DEFAULT_SYSTEMS = ["diag5", "cyclic6", "katsura5"]


def build(name):
    """Return (system, expected_count) for a registered name."""
    if name not in SYSTEMS:
        raise KeyError(f"unknown system {name!r}; known: {', '.join(sorted(SYSTEMS))}")
    builder, expected = SYSTEMS[name]
    return builder(), expected
