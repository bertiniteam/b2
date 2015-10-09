//This file is part of Bertini 2.0.
//
//powerseries_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//powerseries_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with euler_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  powerseries_test.cpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015



#include <boost/test/unit_test.hpp>
#include "tracking/powerseries.hpp"

#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"
#include "tracking/correct.hpp"
#include <deque>

using System = bertini::System;
using Variable = bertini::node::Variable;
using Endgame = bertini::tracking::endgame;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = boost::multiprecision::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(powerseries_endgame_basics)

BOOST_AUTO_TEST_CASE(Hermite_Test_Case)
{
	std::deque< Vec<dbl> > times;

	Vec<dbl> times(1);
	times << dbl(-1);
	

	Vec<dbl> samples(1);
	current_space << ;


	BOOST_CHECK("implemented case where newton step requests higher precision due to AMP criterion B"=="true");

}


BOOST_AUTO_TEST_SUITE_END()
