"""Every conversion out of a multiprecision number, for all four of the types (#389).

Conversion is where a multiprecision value meets code that only knows about floats and
ints -- plotting libraries, ``pandas``, ``json``, anything that calls ``float()``.  Each
conversion is pinned here: what it returns, where it truncates, and where it refuses.

Three of these used to be worse than wrong.  ``int()`` on ``real_mp`` or ``complex_mp``
and ``byteswap()`` on either segfaulted the interpreter, because both types inherit from
``numpy.generic`` once their dtypes are registered, and the inherited implementations
reach for machinery that a user dtype does not have.  ``bool()`` on ``int_mp`` and
``rational_mp`` returned ``True`` for zero, because neither had a ``__bool__`` at all and
python's "every object is true" default applied.  The last test in this file guards the
whole family by exercising every inherited member in a subprocess.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import real_mp, complex_mp, int_mp, rational_mp


BIG = '1' + '0' * 400     # far past what a double can hold


# ----- real_mp ---------------------------------------------------------------------

def test_float_of_real_mp_is_the_nearest_double():
    assert float(real_mp('3.5')) == 3.5
    assert abs(float(real_mp(1) / real_mp(3)) - 1.0 / 3.0) < 1e-16


def test_float_of_real_mp_overflows_to_infinity():
    # float-like semantics, deliberately unlike the exact types below: python's own
    # float(Decimal('1e400')) is inf too, while float(10**400) raises OverflowError
    assert float(real_mp('1e400')) == float('inf')


def test_int_of_real_mp_truncates_toward_zero():
    assert int(real_mp('3.9')) == 3
    assert int(real_mp('-3.9')) == -3
    assert int(real_mp('-0.5')) == 0


def test_int_of_real_mp_keeps_every_digit_before_the_point():
    assert int(real_mp('1e40')) == 10 ** 40
    assert int(real_mp(str(2 ** 80))) == 2 ** 80


def test_int_of_real_mp_refuses_nan_and_infinity():
    with pytest.raises(ValueError):
        int(real_mp('nan'))
    with pytest.raises(OverflowError):
        int(real_mp('inf'))
    with pytest.raises(OverflowError):
        int(real_mp('-inf'))


def test_complex_and_bool_of_real_mp():
    assert complex(real_mp('3.5')) == complex(3.5, 0)
    assert bool(real_mp('0.5')) is True
    assert bool(real_mp(0)) is False


# ----- complex_mp ------------------------------------------------------------------

def test_complex_of_complex_mp_keeps_both_parts():
    assert complex(complex_mp(real_mp('3.5'), real_mp('1.25'))) == complex(3.5, 1.25)


def test_complex_mp_refuses_float_and_int_the_way_python_complex_does():
    z = complex_mp(real_mp('3.5'), real_mp('1.25'))

    with pytest.raises(TypeError) as as_float:
        float(z)
    with pytest.raises(TypeError) as as_int:
        int(z)

    # and each says what to reach for instead
    for raised in (as_float, as_int):
        assert '.real' in str(raised.value)


def test_bool_of_complex_mp():
    assert bool(complex_mp(0, 1)) is True
    assert bool(complex_mp(0, 0)) is False


# ----- int_mp ----------------------------------------------------------------------

def test_int_of_int_mp_is_exact_at_any_size():
    assert int(int_mp(7)) == 7
    assert int(int_mp(BIG)) == 10 ** 400


def test_float_of_int_mp_raises_rather_than_returning_infinity():
    # python's own float(10**400) raises OverflowError; an exact integer type matches it
    assert float(int_mp(7)) == 7.0
    with pytest.raises(OverflowError):
        float(int_mp(BIG))


def test_int_mp_serves_as_an_index():
    assert ['a', 'b', 'c'][int_mp(2)] == 'c'
    assert list(range(int_mp(3))) == [0, 1, 2]


def test_complex_and_bool_of_int_mp():
    assert complex(int_mp(7)) == complex(7, 0)
    assert bool(int_mp(7)) is True
    assert bool(int_mp(0)) is False


# ----- rational_mp -----------------------------------------------------------------

def test_float_of_rational_mp_is_the_nearest_double():
    assert abs(float(rational_mp(1, 3)) - 1.0 / 3.0) < 1e-16


def test_float_of_rational_mp_raises_rather_than_returning_infinity():
    with pytest.raises(OverflowError):
        float(rational_mp(int_mp(BIG), int_mp(1)))


def test_int_of_rational_mp_truncates_toward_zero():
    assert int(rational_mp(22, 7)) == 3
    assert int(rational_mp(-22, 7)) == -3


def test_complex_and_bool_of_rational_mp():
    assert complex(rational_mp(1, 2)) == complex(0.5, 0)
    assert bool(rational_mp(1, 2)) is True
    assert bool(rational_mp(0, 1)) is False


# ----- array-level casts -----------------------------------------------------------

def test_real_mp_array_casts():
    xs = np.array([real_mp(1) / real_mp(4), real_mp(3) / real_mp(4)])
    assert xs.dtype == np.dtype(real_mp)

    assert list(xs.astype(float)) == [0.25, 0.75]
    assert list(xs.astype(complex)) == [complex(0.25, 0), complex(0.75, 0)]
    assert all(isinstance(v, real_mp) for v in xs.astype(object))


def test_complex_mp_array_casts():
    zs = np.array([complex_mp(1, 2), complex_mp(3, 4)])
    assert zs.dtype == np.dtype(complex_mp)

    assert list(zs.astype(complex)) == [complex(1, 2), complex(3, 4)]
    assert all(isinstance(v, complex_mp) for v in zs.astype(object))


def test_object_arrays_of_the_exact_types_cast_to_float():
    # int_mp and rational_mp have no registered dtype, so they land in object arrays;
    # the cast runs element by element through __float__
    assert list(np.array([int_mp(1), int_mp(2)]).astype(float)) == [1.0, 2.0]
    assert list(np.array([rational_mp(1, 2)]).astype(float)) == [0.5]


def test_there_is_no_automatic_promotion_to_float64():
    """Casting is explicit, always: an automatic one would silently drop digits."""
    real = np.dtype(real_mp)

    assert not np.can_cast(real, float, 'safe')
    assert not np.can_cast(real, float, 'same_kind')
    assert not np.can_cast(float, real, 'safe')
    assert np.can_cast(real, float, 'unsafe')       # what .astype(float) uses


# ----- the inherited numpy.generic members -----------------------------------------

def test_byteswap_refuses_on_a_scalar():
    # the value is a handle to digits held elsewhere; there is no byte order to swap
    with pytest.raises(TypeError):
        real_mp(1).byteswap()
    with pytest.raises(TypeError):
        complex_mp(1, 2).byteswap()


def test_byteswap_leaves_an_array_alone():
    xs = np.array([real_mp(1), real_mp(2)])
    swapped = xs.byteswap()
    assert [str(v) for v in swapped] == [str(v) for v in xs]


def test_no_inherited_numpy_member_crashes_the_interpreter():
    """Call every inherited ``numpy.generic`` member; a crash would kill the child only.

    ``real_mp`` and ``complex_mp`` become subclasses of ``numpy.generic`` when their
    dtypes are registered, which hands them about a hundred methods implemented against
    numpy's own scalars.  Several of those reach for cast functions or data layouts a
    user dtype does not have.  Raising is fine -- dying is not, and a segfault inside
    pytest would take every other test down with it, so this runs out of process.
    """
    child = textwrap.dedent("""
        import numpy as np
        import bertini
        from bertini.multiprec import real_mp, complex_mp

        bertini.default_precision(50)

        for value in (real_mp('3.5'), complex_mp(real_mp('3.5'), real_mp('1.25'))):
            own = set(type(value).__dict__)
            for name in sorted(n for n in dir(np.generic) if n not in own):
                try:
                    member = getattr(value, name)
                    if callable(member):
                        member()
                except Exception:
                    pass
        print('SWEEP-OK')
    """)

    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    finished = subprocess.run([sys.executable, '-c', child],
                              capture_output=True, text=True, env=env, timeout=300)

    assert finished.returncode == 0, (
        "a numpy.generic member killed the interpreter (exit %s)\n%s"
        % (finished.returncode, finished.stderr[-2000:]))
    assert 'SWEEP-OK' in finished.stdout
