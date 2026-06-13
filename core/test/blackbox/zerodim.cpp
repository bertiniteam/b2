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


// The CLI infers the start system from the variable-group structure (classic
// Bertini behavior), rather than reading a setting: a single affine variable
// group -> total degree; multiple groups or any homogeneous group -> mhom.
BOOST_AUTO_TEST_SUITE(start_system_inference)

using namespace bertini;
using Var = std::shared_ptr<bertini::node::Variable>;

BOOST_AUTO_TEST_CASE(single_affine_group_infers_total_degree)
{
	Var x = node::Variable::Make("x");
	Var y = node::Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y - 1);
	sys.AddFunction(x + y);

	BOOST_CHECK(blackbox::InferStartType(sys) == blackbox::type::Start::TotalDegree);
}

BOOST_AUTO_TEST_CASE(two_affine_groups_infer_mhom)
{
	Var x = node::Variable::Make("x");
	Var y = node::Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x - y);

	BOOST_CHECK(blackbox::InferStartType(sys) == blackbox::type::Start::MHom);
}

BOOST_AUTO_TEST_CASE(homogeneous_group_infers_mhom)
{
	// a single homogeneous variable group is projective, not total-degree
	Var x = node::Variable::Make("x");
	Var y = node::Variable::Make("y");
	System sys;
	sys.AddHomVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y);

	BOOST_CHECK(blackbox::InferStartType(sys) == blackbox::type::Start::MHom);
}

BOOST_AUTO_TEST_CASE(griewank_osborn_single_group_infers_total_degree)
{
	auto sys = system::Precon::GriewankOsborn();
	BOOST_CHECK(blackbox::InferStartType(sys) == blackbox::type::Start::TotalDegree);
}

// DISABLED pending the block-composed MHom start system (plan
// typed-wondering-sutherland, tasks #4-#8).  Today MakeZeroDim for an MHom
// target throws "unknown visitor: LinearProduct" because the start system is
// built from LinearProduct nodes the SLP compiler can't compile.  Once the MHom
// start is rebuilt on ProductsOfLinearsBlock, re-enable this and extend it from
// "constructs" to "actually solves".
BOOST_AUTO_TEST_CASE(inferred_mhom_builds_an_mhomogeneous_zerodim,
                     * boost::unit_test::disabled())
{
	// end to end: a two-variable-group system, inferred to MHom, must build a
	// ZeroDim backed by the MHomogeneous start system (not total degree).
	using namespace bertini::tracking;
	using namespace bertini::endgame;

	Var x = node::Variable::Make("x");
	Var y = node::Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x - y);

	blackbox::ZeroDimRT rt;
	rt.start = blackbox::InferStartType(sys);
	rt.tracker = blackbox::type::Tracker::Adaptive;
	rt.endgame = blackbox::type::Endgame::Cauchy;
	BOOST_CHECK(rt.start == blackbox::type::Start::MHom);

	using MHomZD = algorithm::ZeroDim<AMPTracker, typename EndgameSelector<AMPTracker>::Cauchy, System, start_system::MHomogeneous>;
	using TotDegZD = algorithm::ZeroDim<AMPTracker, typename EndgameSelector<AMPTracker>::Cauchy, System, start_system::TotalDegree>;

	auto zd = blackbox::MakeZeroDim(rt, sys);
	BOOST_CHECK(dynamic_cast<MHomZD*>(zd.get()) != nullptr);
	BOOST_CHECK(dynamic_cast<TotDegZD*>(zd.get()) == nullptr);
}

BOOST_AUTO_TEST_SUITE_END() // end the start_system_inference sub-suite

BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite
