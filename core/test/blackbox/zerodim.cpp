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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame

#include <boost/test/unit_test.hpp>

#include "bertini2/blackbox/switches_zerodim.hpp"
#include "bertini2/system/precon.hpp"


BOOST_AUTO_TEST_SUITE(blackbox_test)

BOOST_AUTO_TEST_SUITE(zero_dim)

using namespace bertini;

BOOST_AUTO_TEST_CASE(make_zero_dim_defaults)
{
	auto sys = system::Precon::GriewankOsborn();
	blackbox::ZeroDimRT my_runtime_type_options; // make defaults

	auto zd_ptr = blackbox::MakeZeroDim(my_runtime_type_options, sys);
}



BOOST_AUTO_TEST_CASE(make_zero_dim_nondefaults)
{
	auto sys = system::Precon::GriewankOsborn();
	blackbox::ZeroDimRT my_runtime_type_options; // make defaults
	my_runtime_type_options.start = bertini::blackbox::type::Start::MHom;
	my_runtime_type_options.tracker = bertini::blackbox::type::Tracker::FixedDouble;
	my_runtime_type_options.endgame = bertini::blackbox::type::Endgame::Cauchy;
	auto zd_ptr = blackbox::MakeZeroDim(my_runtime_type_options, sys);


	// zd_ptr->DefaultSetup();
}


BOOST_AUTO_TEST_CASE(make_zero_dim_honors_endgame_choice)
{
	// regression: ZeroDimSpecifyShouldClone hardcoded the Cauchy endgame,
	// ignoring its EndgameType parameter -- selecting PowerSeries silently
	// built a Cauchy ZeroDim.  pin the dynamic types.
	using namespace bertini::tracking;
	using namespace bertini::endgame;
	auto sys = system::Precon::GriewankOsborn();

	blackbox::ZeroDimRT rt;
	rt.tracker = blackbox::type::Tracker::Adaptive;

	using PSEGZD   = algorithm::ZeroDim<AMPTracker, typename EndgameSelector<AMPTracker>::PSEG,   System, start_system::TotalDegree>;
	using CauchyZD = algorithm::ZeroDim<AMPTracker, typename EndgameSelector<AMPTracker>::Cauchy, System, start_system::TotalDegree>;

	rt.endgame = blackbox::type::Endgame::PowerSeries;
	auto zd_pseg = blackbox::MakeZeroDim(rt, sys);
	BOOST_CHECK(dynamic_cast<PSEGZD*>(zd_pseg.get()) != nullptr);
	BOOST_CHECK(dynamic_cast<CauchyZD*>(zd_pseg.get()) == nullptr);

	rt.endgame = blackbox::type::Endgame::Cauchy;
	auto zd_cauchy = blackbox::MakeZeroDim(rt, sys);
	BOOST_CHECK(dynamic_cast<CauchyZD*>(zd_cauchy.get()) != nullptr);
	BOOST_CHECK(dynamic_cast<PSEGZD*>(zd_cauchy.get()) == nullptr);
}

BOOST_AUTO_TEST_SUITE_END() // end the zerodim sub-suite

BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite
