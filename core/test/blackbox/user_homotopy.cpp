//This file is part of Bertini 2.
//
//test/blackbox/zerodim.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/blackbox/zerodim.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/blackbox/zerodim.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame

#include <boost/test/unit_test.hpp>

#include "bertini2/blackbox/switches_zerodim.hpp"
#include "bertini2/system/precon.hpp"

#include "bertini2/blackbox/user_homotopy.hpp"
BOOST_AUTO_TEST_SUITE(blackbox_test)

BOOST_AUTO_TEST_SUITE(user_hom)

using namespace bertini;




BOOST_AUTO_TEST_CASE(make_user_hom_dim_nondefaults)
{
	auto sys1 = system::Precon::GriewankOsborn();
	SampCont<mpfr> start_solns;

	auto sys2 = start_system::User(sys1, start_solns);
	auto sys3 = system::Precon::GriewankOsborn();

	blackbox::ZeroDimRT my_runtime_type_options; // make defaults

	my_runtime_type_options.start = bertini::blackbox::type::Start::User;
	my_runtime_type_options.tracker = bertini::blackbox::type::Tracker::FixedDouble;
	my_runtime_type_options.endgame = bertini::blackbox::type::Endgame::Cauchy;
	auto zd_ptr = blackbox::MakeUserHom(my_runtime_type_options, sys1, sys2, sys3);

	
	// zd_ptr->DefaultSetup();
}

BOOST_AUTO_TEST_SUITE_END() // end the user_hom sub-suite

BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite
