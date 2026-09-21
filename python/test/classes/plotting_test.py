"""Plotting multiprecision numbers with matplotlib, and the one thing it refuses (#389).

``real_mp`` data goes straight onto a pyplot axis.  The calls that compute before they
draw -- ``hist``, ``bar`` with a plain float width -- need a units converter, because
there is no common dtype between ``real_mp`` and ``float64`` and deliberately never will
be.  :mod:`bertini._matplotlib_bridge` registers one at import, so this file is mostly a
list of calls that must keep working.

The refusals matter more than the successes.  numpy's ``.real`` / ``.imag`` on an array
of a user complex dtype return the array itself and zeros, so ``scatter(z.real, z.imag)``
would put every point on the x axis with no error -- which is why complex data is turned
away at the axis instead.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import real_mp, complex_mp

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402  (after the backend is chosen)


@pytest.fixture
def axes():
    figure, axis = plt.subplots()
    yield axis
    plt.close(figure)


@pytest.fixture
def xs():
    return np.array([real_mp(k) / real_mp(7) for k in range(1, 6)])


@pytest.fixture
def zs():
    return np.array([complex_mp(real_mp(k) / real_mp(7), real_mp(k))
                     for k in range(1, 6)])


# ----- what must simply work -------------------------------------------------------

def test_the_ordinary_calls_take_mp_arrays(axes, xs):
    ys = xs * xs

    axes.plot(xs, ys)
    axes.scatter(xs, ys)
    axes.step(xs, ys)
    axes.loglog(xs, ys)
    axes.fill_between(xs, ys, xs)
    axes.errorbar(xs, ys, yerr=xs)
    axes.stem(xs, ys)
    axes.axhline(xs[0])
    axes.axvline(xs[0])
    axes.text(xs[0], ys[0], 'a label')
    axes.set_xlim(xs[0], xs[-1])

    axes.figure.savefig(os.devnull, format='png')


def test_the_ordinary_calls_take_lists_too(axes, xs):
    values = list(xs)
    axes.plot(values, values)
    axes.scatter(values, values)


def test_histogram_and_bar_need_the_converter(axes, xs):
    """These compute against plain floats internally, so they fail without one."""
    axes.hist(xs)
    axes.hist(list(xs))
    axes.bar(xs, xs, width=0.01)


def test_values_arrive_where_they_belong(axes, xs):
    drawn = axes.scatter(xs, xs * xs).get_offsets()

    assert [round(x, 12) for x in drawn[:, 0]] == [round(float(v), 12) for v in xs]
    assert [round(y, 12) for y in drawn[:, 1]] == [round(float(v * v), 12) for v in xs]


def test_colour_can_be_mapped_from_mp_values(axes, xs):
    axes.scatter(xs, xs, c=xs)


def test_the_complex_plane_with_the_parts_taken_explicitly(axes, zs):
    drawn = axes.scatter(pb.real(zs), pb.imag(zs)).get_offsets()

    assert [round(y, 12) for y in drawn[:, 1]] == [round(float(z.imag), 12) for z in zs]


def test_a_solution_plots_in_the_complex_plane(axes):
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group([x, y])
    system.add_function(x * x + y * y - 1)
    system.add_function(x - y)

    solution = pb.solve(system).solutions[0]

    # a Solution overrides .real/.imag at the subclass level and is simply correct
    drawn = axes.scatter(solution.real, solution.imag).get_offsets()
    assert [round(y, 12) for y in drawn[:, 1]] == [round(float(v.imag), 12)
                                                   for v in np.asarray(solution)]


# ----- what must be refused --------------------------------------------------------

def test_complex_data_is_refused_at_an_axis(axes, zs):
    with pytest.raises(TypeError) as raised:
        axes.plot(zs)

    assert 'bertini.real' in str(raised.value)


def test_the_silent_trap_is_refused(axes, zs):
    """``.real``/``.imag`` on a plain mp-complex array are the array and zeros."""
    with pytest.raises(TypeError):
        axes.scatter(zs.real, zs.imag)


def test_complex_data_is_refused_even_after_real_data_set_the_axis(axes, xs, zs):
    axes.plot(xs, xs)     # the axis now holds the real converter

    with pytest.raises(TypeError):
        axes.scatter(zs.real, zs.imag)


def test_complex_data_is_refused_by_the_computing_calls_too(axes, zs):
    with pytest.raises(TypeError):
        axes.hist(zs)


# ----- registration ----------------------------------------------------------------

def _run(child):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p),
               MPLBACKEND='Agg')
    return subprocess.run([sys.executable, '-c', textwrap.dedent(child)],
                          capture_output=True, text=True, env=env, timeout=300)


@pytest.mark.parametrize('bertini_first', [True, False])
def test_the_converter_registers_whichever_import_comes_first(bertini_first):
    imports = ("import bertini\n        import matplotlib.pyplot as plt"
               if bertini_first else
               "import matplotlib.pyplot as plt\n        import bertini")

    finished = _run("""
        %s
        import numpy as np
        from bertini.multiprec import real_mp

        figure, axis = plt.subplots()
        axis.hist(np.array([real_mp(k) / real_mp(7) for k in range(1, 6)]))
        print('HIST-OK')
    """ % imports)

    assert finished.returncode == 0, finished.stderr[-2000:]
    assert 'HIST-OK' in finished.stdout


def test_importing_bertini_does_not_require_matplotlib():
    finished = _run("""
        import sys

        class Blocker:
            def find_spec(self, fullname, path=None, target=None):
                if fullname == 'matplotlib' or fullname.startswith('matplotlib.'):
                    raise ImportError('no matplotlib here')
                return None

        sys.meta_path.insert(0, Blocker())

        import bertini
        from bertini.multiprec import real_mp
        assert float(real_mp('0.5')) == 0.5
        assert 'matplotlib' not in sys.modules
        print('NO-MATPLOTLIB-OK')
    """)

    assert finished.returncode == 0, finished.stderr[-2000:]
    assert 'NO-MATPLOTLIB-OK' in finished.stdout


# ----- pandas, which draws through matplotlib --------------------------------------

def test_pandas_scatter_draws_mp_columns():
    pandas = pytest.importorskip("pandas")

    frame = pandas.DataFrame({'x': [real_mp(k) / real_mp(7) for k in range(1, 6)],
                              'y': [real_mp(k) for k in range(1, 6)]})
    axis = frame.plot.scatter(x='x', y='y')
    plt.close(axis.figure)


def test_pandas_screens_by_dtype_before_matplotlib_sees_the_data():
    """``.plot()`` and ``.plot.hist()`` select numeric columns themselves first.

    mp values live in an object column, so pandas drops them before any converter can
    run.  Casting the column is the way through, and it is exact up to the double.
    """
    pandas = pytest.importorskip("pandas")

    frame = pandas.DataFrame({'x': [real_mp(k) / real_mp(7) for k in range(1, 6)],
                              'y': [real_mp(k) for k in range(1, 6)]})

    with pytest.raises(TypeError):
        frame.plot(x='x', y='y')

    axis = frame.astype(float).plot(x='x', y='y')
    plt.close(axis.figure)
