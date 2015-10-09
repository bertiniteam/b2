//This file is part of Bertini 2.0.
//
//tracker_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracker_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracker_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015

//tracker_test.cpp
//


#include <boost/test/unit_test.hpp>

#include "tracking/tracker.hpp"

using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(tracker_basics)



BOOST_AUTO_TEST_CASE(tracker_track_linear)
{
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");

	System sys;

	sys.AddFunction(x-t);
	sys.AddPathVariable(t);


	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	config::Stepping stepping_preferences;
	config::Newton newton_preferences;
	config::AdaptiveMultiplePrecisionConfig AMP;

	tracker.Setup(config::Predictor::Euler,
					mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0.0");
	
	Vec<mpfr> x_start(1);
	x_start << mpfr("1.0");

	Vec<mpfr> x_end(1);

	tracker.TrackPath(x_end,
	                  t_start, t_end, x_start);

	BOOST_CHECK_EQUAL(x_end.size(),1);

	BOOST_CHECK(abs(x_end(1)-mpfr("1.0")) < 1e-5);

}



BOOST_AUTO_TEST_SUITE_END()




