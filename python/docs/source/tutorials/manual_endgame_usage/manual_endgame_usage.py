"""Tutorial: Using an endgame to compute singular endpoints (manual endgame usage).

Assembles every ``.. testcode::`` fragment from index.rst into one runnable program.
Run:  python manual_endgame_usage.py
"""

import bertini


def form_system():
    """Build the Griewank-Osborne system from scratch."""
    gw = bertini.System()

    x = bertini.Variable("x")
    y = bertini.Variable("y")

    vg = bertini.VariableGroup()
    vg.append(x)
    vg.append(y)
    gw.add_variable_group(vg)

    gw.add_function(bertini.multiprec.rational_mp(29, 16) * x**3 - 2 * x * y)
    gw.add_function(y - x**2)

    return gw


def form_homotopy(gw):
    """Make the total degree start system and couple it via the gamma trick."""
    t = bertini.Variable('t')
    td = bertini.system.start_system.TotalDegreeLinearProduct(gw)
    gamma = bertini.symbolics.Rational.rand()
    hom = (1 - t) * gw + t * gamma * td
    hom.add_path_variable(t)

    return hom, td


def track_to_boundary(gw, hom, td):
    """Track every start point to the endgame boundary."""
    tr = bertini.AMPTracker(hom)

    start_time = bertini.multiprec.complex_mp("1")
    eg_boundary = bertini.multiprec.complex_mp("0.1")

    midpath_points = [None] * td.num_start_points()
    for ii in range(td.num_start_points()):
        midpath_points[ii] = bertini.multiprec.Vector(gw.num_variables())   # result must be pre-sized
        code = tr.track_path(result=midpath_points[ii], start_time=start_time, end_time=eg_boundary, start_point=td.start_point_mp(ii))
        assert code == bertini.SuccessCode.Success                 # every path reaches the boundary

    return tr, eg_boundary, midpath_points


def use_endgame(tr, eg_boundary, td, midpath_points):
    """Run the adaptive Cauchy endgame from the boundary down to t = 0."""
    eg = bertini.endgame.AMPCauchyEndgame(tr, eg_boundary)

    # make an observer to be able to see what's going on inside
    ob = bertini.endgame.observers.amp_cauchy.GoryDetailLogger()

    eg.add_observer(ob)

    # nothing has run yet, so things are empty and default
    assert eg.cycle_number() == 0
    assert len(eg.final_approximation()) == 0    # nothing computed yet

    final_points = []
    codes = []
    for ii in range(td.num_start_points()):
        codes.append(eg.run(midpath_points[ii]))   # refine from the boundary down to t = 0
        final_points.append(eg.final_approximation())

    # exactly three of the six paths converge to the triple point at the origin
    origin_hits = sum(1 for fa in final_points
                      if len(fa) and max(abs(complex(v)) for v in fa) < 1e-6)
    print('paths landing on the triple point:', origin_hits)

    return final_points, origin_hits


def main():
    gw = form_system()
    hom, td = form_homotopy(gw)
    tr, eg_boundary, midpath_points = track_to_boundary(gw, hom, td)
    use_endgame(tr, eg_boundary, td, midpath_points)


if __name__ == '__main__':
    main()
