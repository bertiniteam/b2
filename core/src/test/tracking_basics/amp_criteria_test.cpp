//This file is part of Bertini 2.0.
//
//amp_criteria_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_criteria_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_criteria_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  amp_criteria_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015



#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"

#include "tracking/amp_criteria.hpp"


using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(amp_criteria_tracking_basics)


BOOST_AUTO_TEST_CASE(AMP_criteriaA_double)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaA_mp)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}


BOOST_AUTO_TEST_CASE(AMP_criteriaB_double)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaB_mp)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}

BOOST_AUTO_TEST_CASE(AMP_criteriaC_double)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaC_mp)
{
	BOOST_CHECK_EQUAL("test implemented","true");
}




BOOST_AUTO_TEST_SUITE_END()

