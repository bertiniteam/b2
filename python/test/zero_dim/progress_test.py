"""Reporting how far a solve has got (#359).

The display is on by default, so most of what matters here is that it stays *invisible*:
a solve inside a script, a test or a pipeline must print nothing, and a solve that finds
its answers in the records must not show a bar that never moves.

The counting is checked separately from the drawing, because the counting is the part
that can be wrong in a way nobody notices.
"""

import io

import pytest

import bertini as pb
import bertini._progress as progress


@pytest.fixture
def ordinary():
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group(pb.VariableGroup([x, y]))
    system.add_function(x ** 3 + y ** 3 - 1)
    system.add_function(x ** 3 - y ** 2 + x)
    return system


@pytest.fixture
def has_a_path_to_infinity():
    """y = x^2 meeting xy = 1: Bezout says four paths, only three endpoints are finite."""
    x, y = pb.Variable('x'), pb.Variable('y')
    system = pb.System()
    system.add_variable_group(pb.VariableGroup([x, y]))
    system.add_function(x * x - y)
    system.add_function(x * y - 1)
    return system


def _watched(system, drawn_to=None, threads=1, **kwargs):
    """Solve while watching, and hand back the report."""
    pb.recording(False)
    solver = pb.ZeroDimSolver(system, mptype='adaptive')
    if threads != 1:
        solver.update(num_threads=threads)

    report = progress.ProgressReport(solver, stream=drawn_to or io.StringIO(), **kwargs)
    solver.add_observer(report)
    solver.solve(show_progress=False)      # the one being tested is ours, attached above
    solver.remove_observer(report)
    return report


def test_every_path_is_counted_once(ordinary):
    report = _watched(ordinary)

    assert report.total == 9                      # Bezout: 3 x 3
    assert report.finished == report.total
    assert report.tallies == {'ok': 9}


def test_every_path_is_counted_once_under_threads(ordinary):
    """Path events arrive from worker threads, out of order."""
    report = _watched(ordinary, threads=4)

    assert report.total == 9
    assert report.finished == 9
    assert sum(report.tallies.values()) == 9


def test_outcomes_are_tallied_as_they_arrive(has_a_path_to_infinity):
    report = _watched(has_a_path_to_infinity)

    assert report.tallies.get('diverged') == 1
    assert report.tallies.get('ok') == 3
    assert report.tally_text() == '3 ok, 1 diverged'


def test_nothing_is_drawn_when_the_stream_is_not_a_terminal(ordinary):
    drawn = io.StringIO()
    report = _watched(ordinary, drawn_to=drawn, delay=0.0)

    assert report.finished == 9          # it counted
    assert drawn.getvalue() == ''        # and drew nothing


def test_forcing_it_draws(ordinary):
    drawn = io.StringIO()
    _watched(ordinary, drawn_to=drawn, delay=0.0, force=True)

    written = drawn.getvalue()
    assert 'paths' in written
    assert '9/9' in written
    assert '9 ok' in written


def test_a_recalled_solve_reports_recall_instead_of_a_bar(ordinary, tmp_path):
    directory = str(tmp_path / 'records')
    pb.recording(True)

    def one_solve(watched):
        pb.random.set_random_seed(7)
        solver = pb.ZeroDimSolver(ordinary, mptype='adaptive')
        solver.record_to(directory)
        report = None
        if watched:
            report = progress.ProgressReport(solver, stream=io.StringIO(),
                                             delay=0.0, force=True)
            solver.add_observer(report)
        solver.solve(show_progress=False)
        if watched:
            solver.remove_observer(report)
        return solver, report

    first, _ = one_solve(watched=False)
    assert first.num_paths_recalled() == 0

    second, report = one_solve(watched=True)
    assert second.num_paths_recalled() == 9, 'the second solve was not the same ask'

    assert report.recalled == 9
    assert report.total == 9
    assert report.finished == 0          # nothing was tracked, so no bar advanced
    assert 'recalled 9 of 9 paths' in report.stream.getvalue()


def test_the_redraw_interval_scales_with_the_path_count():
    """b1 called this printpathprogress: at a million paths, drawing on each one is the problem."""

    class Solver:
        def __init__(self, paths):
            self.paths = paths

        def num_paths(self):
            return self.paths

        def num_paths_recalled(self):
            return 0

    for paths, expected in ((9, 1), (100, 1), (10_000, 50), (1_000_000, 5000)):
        report = progress.ProgressReport(Solver(paths))
        report._start()
        report._size_to_the_work()
        assert report.every == expected, 'at %d paths' % paths
        assert paths // max(report.every, 1) <= progress.TARGET_UPDATES + 1


def test_an_explicit_interval_is_kept():
    report = progress.ProgressReport(None, every=25)
    assert report.every == 25


def test_the_report_never_breaks_a_solve(ordinary):
    """A display must not be the reason a solve fails."""

    class Broken(progress.ProgressReport):
        def _observe(self, event):
            raise RuntimeError('the display is broken')

    pb.recording(False)
    solver = pb.ZeroDimSolver(ordinary, mptype='adaptive')
    broken = Broken(solver, stream=io.StringIO())
    solver.add_observer(broken)
    solver.solve(show_progress=False)
    solver.remove_observer(broken)

    assert len(solver.all_solutions()) == 9


def test_show_progress_false_attaches_nothing(ordinary):
    pb.recording(False)
    solver = pb.ZeroDimSolver(ordinary, mptype='adaptive')

    with progress.watching(solver, show_progress=False) as report:
        assert report is None

    result = solver.solve(show_progress=False)
    assert len(result.answer) == 9


def test_solving_through_the_casual_api_stays_silent(ordinary, capsys):
    """show_progress defaults to True; a test run is not a terminal, so nothing appears."""
    pb.recording(False)
    result = pb.solve(ordinary)

    captured = capsys.readouterr()
    assert len(result.answer) == 9
    assert captured.out == ''
    assert captured.err == ''
