"""Bertini 2 tutorial: Tracking to nonsingular endpoints.

The complete, runnable assembly of the tutorial's testcode fragments.
Run:  python tracking_nonsingular.py
"""

import bertini
import numpy as np


def form_system():
    """Build the target polynomial system sys = {x^2 + y^2 - 1, x + y}."""
    x = bertini.symbolics.Variable("x")  # variable name need not match
    y = bertini.symbolics.Variable("y")

    f = x**2 + y**2 - 1  # ** is exponentiation in Python.
    g = x + y

    sys = bertini.System()
    sys.add_function(f)
    sys.add_function(g)

    # group the variables affinely and stuff them into sys
    grp = bertini.VariableGroup()
    grp.append(x)
    grp.append(y)
    sys.add_variable_group(grp)

    # check the degrees of our functions
    d = sys.degrees()
    assert(d[0] == 2)  # f is degree 2 (highest power in any term is 2)
    assert(d[1] == 1)  # g is degree 1 (highest power in any term is 1)

    return sys, grp, f, g


def explore_nonalgebraic(sys, grp, f, g):
    """Aside: adding non-polynomial functions, then rebuild an algebraic sys."""
    x = grp[0]  # recover x from the variable group for the aside
    sys.add_function(x**-1)  # happily accepts a non-polynomial function.
    sys.add_function(bertini.symbolics.sin(x))
    d = sys.degrees()
    assert(d[2] == -1)  # unsurprising, but actually a coincidence
    assert(d[3] == -1)  # also -1.  anything non-polynomial is a negative number.
    # sin has no well-defined degree
    # bertini uses negative degree to indicate non-polynomial

    # correcting our system -- a return to algebraicness.
    del sys  # we mal-formed our system above, so we start over
    sys = bertini.System()
    sys.add_variable_group(grp)
    sys.add_function(f)
    sys.add_function(g)
    return sys


def form_start_system(sys):
    """Make a Total Degree start system mirroring sys, and sample start points."""
    td = bertini.system.start_system.TotalDegreeLinearProduct(sys)

    # generate the 1th (0-based offsets in python) start point
    sp_d = td.start_point_d(1)  # at double precision

    sp_mp = td.start_point_mp(1)  # at current default multiple precision
    assert(bertini.default_precision() == sp_mp[1].precision)

    return td


def form_homotopy(sys, td):
    """Glue sys and td into a homotopy in the path variable t."""
    t = bertini.Variable("t")     # make a path variable
    homotopy = (1 - t) * sys + t * td   # glue
    homotopy.add_path_variable(t)   # indicate the path var
    return homotopy


def track_single_path(homotopy, td):
    """Configure an AMP tracker and track a single path, with logging demo."""
    tr = bertini.AMPTracker(homotopy)
    tr.tracking_tolerance(1e-5)  # track the path to 5 digits or so

    # adjust some stepping settings
    stepping = bertini.tracking.SteppingConfig()
    stepping.max_step_size = bertini.multiprec.real_mp(1) / bertini.multiprec.real_mp(13)

    # then, set the config into the tracker.
    tr.set_stepping(stepping)

    result = np.zeros((2,), dtype=bertini.multiprec.complex_mp)
    tr.track_path(result, bertini.multiprec.complex_mp(1),
                  bertini.multiprec.complex_mp(0), td.start_point_mp(0))

    # make an observer to log details, and attach it
    g = bertini.tracking.observers.amp.GoryDetailLogger()
    tr.add_observer(g)

    # re-run: prints a ton of detail to the screen
    tr.track_path(result, bertini.multiprec.complex_mp(1),
                  bertini.multiprec.complex_mp(0), td.start_point_mp(0))

    # turn logging back off
    tr.remove_observer(g)

    return tr, result


def track_all_paths():
    """The 'grab the whole thing' block: track every start path."""
    x = bertini.symbolics.Variable("x")  # variable name need not match
    y = bertini.symbolics.Variable("y")
    f = x**2 + y**2 - 1
    g = x + y

    sys = bertini.System()
    sys.add_function(f)
    sys.add_function(g)

    grp = bertini.VariableGroup()
    grp.append(x)
    grp.append(y)
    sys.add_variable_group(grp)

    td = bertini.system.start_system.TotalDegreeLinearProduct(sys)

    t = bertini.Variable("t")
    homotopy = (1 - t) * sys + t * td
    homotopy.add_path_variable(t)

    tr = bertini.AMPTracker(homotopy)

    # commented out for screen-saving.
    # g = bertini.tracking.observers.amp.GoryDetailLogger()
    # tr.add_observer(g)
    # one could also bertini.logging.init() and set a file name,
    # so it gets piped there instead of wherever Boost.Log goes by default.

    tr.tracking_tolerance(1e-5)  # track the path to 5 digits or so
    tr.infinite_truncation_tolerance(1e5)
    tr.predictor(bertini.Predictor.RK4)
    stepping = bertini.tracking.SteppingConfig()
    stepping.max_step_size = bertini.multiprec.real_mp(1) / bertini.multiprec.real_mp(13)

    # set the config into the tracker
    tr.set_stepping(stepping)

    results = []  # make an empty list into which to put the results
    expected_code = bertini.SuccessCode.Success
    codes = []
    for ii in range(td.num_start_points()):
        results.append(np.zeros((2,), dtype=bertini.multiprec.complex_mp))
        codes.append(tr.track_path(result=results[-1],
                                   start_time=bertini.multiprec.complex_mp(1),
                                   end_time=bertini.multiprec.complex_mp(0),
                                   start_point=td.start_point_mp(ii)))

    # tr.remove_observer(g)

    # the tracked endpoints are now in the list ``results``
    print(codes == [expected_code] * 2)


def main():
    sys, grp, f, g = form_system()
    sys = explore_nonalgebraic(sys, grp, f, g)
    td = form_start_system(sys)
    homotopy = form_homotopy(sys, td)
    track_single_path(homotopy, td)
    track_all_paths()


if __name__ == '__main__':
    main()
