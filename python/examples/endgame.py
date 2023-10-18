import pybertini as pb

from pybertini import Variable, VariableGroup, System
import pybertini.system.start_system as ss

from pybertini.function_tree.symbol import Rational

from pybertini.endgame import *
from pybertini.endgame.config import *

from pybertini.tracking import *
from pybertini.tracking.config import *

from pybertini.multiprec import Float as mpfr_float
from pybertini.multiprec import Complex as mpfr_complex

import numpy as np



ambient_precision = 50;



x = Variable("x");
y = Variable("y");
t = Variable("t");
#
sys = System();
#
var_grp = VariableGroup();
var_grp.append(x);
var_grp.append(y);

sys.add_variable_group(var_grp);

sys.add_function((x-1)**3)
sys.add_function((y-1)**2)

sys.homogenize();
sys.auto_patch();

assert (sys.is_patched() == 1)
assert (sys.is_homogeneous() == 1)

td = ss.TotalDegree(sys);

assert (td.is_patched() == 1)
assert (td.is_homogeneous() == 1)

gamma = Rational.rand();


final_system = (1-t)*sys + gamma*t*td;
final_system.add_path_variable(t);

print(final_system)
prec_config = AMPConfig(final_system);

stepping_pref = SteppingConfig();
newton_pref = NewtonConfig();



tracker = AMPTracker(final_system);

tracker.setup(Predictor.RK4, 1e-5, 1e5, stepping_pref, newton_pref);
tracker.precision_setup(prec_config);

num_paths_to_track = td.num_start_points();
n = int(str(num_paths_to_track));

t_start = mpfr_complex(1);
t_endgame_boundary = mpfr_complex("0.1");
t_final = mpfr_complex(0);


print('tracking to the endgame boundary')


bdry_points = [np.empty(dtype=mpfr_complex, shape=(3,)) for i in range(n)]

for i in range(n):

    pb.default_precision(ambient_precision);
    final_system.precision(ambient_precision);

    print('here')

    td.precision(ambient_precision)
    start_point = td.start_point_mp(i);

    print('there')
    print(pb.multiprec.precision(start_point))


    bdry_pt = np.zeros(dtype=mpfr_complex, shape=(3));
    track_success_code = tracker.track_path(bdry_pt,t_start, t_endgame_boundary, start_point);
    print(bdry_pt.flags)
    bdry_points[i] = bdry_pt;

    assert (track_success_code == SuccessCode.Success)



tracker.setup(Predictor.HeunEuler, 1e-6, 1e5, stepping_pref, newton_pref);
my_endgame = AMPCauchyEG(tracker);


print('running the endgame')

final_homogenized_solutions = [np.empty(dtype=mpfr_complex, shape=(3,)) for i in range(n)]
for i in range(n):
    p = bdry_points[i]


    print('moving to precision {} to match precision of boundary point'.format(pb.multiprec.precision(p)))


    pb.default_precision(p[0].precision());
    final_system.precision(p[0].precision());

    print(p.flags)

    q = np.zeros(dtype=mpfr_complex, shape=(3));
    print(p.flags)
    track_success_code = my_endgame.run(mpfr_complex(t_endgame_boundary),q);
    print('qwfp')
    

    final_homogenized_solutions[i] = my_endgame.final_approximation();
    print(final_system.dehomogenize_point(final_homogenized_solutions[i]));
    assert (track_success_code == SuccessCode.Success)