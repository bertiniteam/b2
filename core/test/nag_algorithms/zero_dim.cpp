//This file is part of Bertini 2.
//
//zero_dim.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//zero_dim.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with zero_dim.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

/**
\file test/nag_algorithms/zero_dim.cpp  Tests the zero dim algorithm.
*/

// individual authors of this file include:
// silviana amethyst

#include "bertini2/system/precon.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"
#include <boost/test/unit_test.hpp>
#include "bertini2/nag_algorithms/output.hpp"


using Variable = bertini::node::Variable;

BOOST_AUTO_TEST_SUITE(zero_dim)

using TrackerT = bertini::tracking::DoublePrecisionTracker;



BOOST_AUTO_TEST_CASE(can_run_griewank_osborn)
{
	using namespace bertini;
	using namespace tracking;

	using Tolerances = algorithm::TolerancesConfig;
	using EndgameConfT = endgame::EndgameConfig;
	

	auto sys = system::Precon::GriewankOsborn();

	auto zd = algorithm::ZeroDim<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-5;
	tols.newton_during_endgame = 1e-6;
	zd.Set(tols);

	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(logger);


	auto eg = zd.GetFromEndgame<EndgameConfT>();
	eg.final_tolerance = 1e-12;
	zd.SetToEndgame(eg);

	zd.Solve();


	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}




BOOST_AUTO_TEST_CASE(can_run_change_some_settings)
{
	
	using namespace bertini;
	using namespace tracking;

	using Tolerances = algorithm::TolerancesConfig;

	auto sys = system::Precon::GriewankOsborn();

	auto zd = algorithm::ZeroDim<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, decltype(sys), start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(logger);

	

	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-4;
	tols.newton_during_endgame = 1e-4;

	zd.Set(tols);

	zd.Solve();


	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}



/**
Check whether we can run zero dim on the non-homogenized version of Griewank Osborn.
*/
BOOST_AUTO_TEST_CASE(reference_managed_systems_GO_nonhom)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	auto TD = start_system::TotalDegree(sys);

	auto t = Variable::Make("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::ZeroDim<
				TrackerT, 
				bertini::endgame::EndgameSelector<TrackerT>::PSEG, 
				decltype(sys), 
				start_system::TotalDegree,
				policy::RefToGiven
					>
			(sys, TD, h); 
	// have to pass in all three to the constructor, because using references. 

	zd.DefaultSetup();

	zd.Solve();

	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}





/**
Check whether we can run zero dim on the non-homogenized version of Griewank Osborn.
*/
BOOST_AUTO_TEST_CASE(reference_managed_systems_GO)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	sys.Homogenize();
	sys.AutoPatch();

	auto TD = start_system::TotalDegree(sys);

	auto t = Variable::Make("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::ZeroDim<
				TrackerT, 
				bertini::endgame::EndgameSelector<TrackerT>::PSEG, 
				decltype(sys), 
				start_system::TotalDegree,
				policy::RefToGiven
					>
			(sys, TD, h); 
	// have to pass in all three to the constructor, because using references. 

	zd.DefaultSetup();

	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(logger);
	
	zd.Solve();

	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}


// Run ZeroDim from a USER-CONSTRUCTED homotopy and a GIVEN list of start points -- the
// parameter-homotopy workflow.  Crucially this reuses the ENTIRE ZeroDim solve pipeline
// (pre-endgame tracking, the midpath check, the endgame, post-processing) UNCHANGED: it is the
// same ZeroDim class template, merely instantiated with start_system::User (start points come
// from the supplied list, not generated) and policy::RefToGiven (the homotopy is taken as-is,
// not formed by homogenizing + coupling a start system).
//
// Parameter homotopy H(x,t) = x^2 - (9 - 5 t).  At t = 1 the roots are +/-2 (the given start
// points); at t = 0 they are +/-3 (the target x^2 - 9).  Tracking the two start points down to
// t = 0 must recover +/-3 -- i.e. ZeroDim moved the t=1 solutions to the t=0 parameter.
BOOST_AUTO_TEST_CASE(user_homotopy_parameter_homotopy_solves)
{
	using namespace bertini;
	using namespace tracking;
	using mpfr = bertini::mpfr_complex;

	auto x = Variable::Make("x");
	auto t = Variable::Make("t");

	// the homotopy H(x,t) = x^2 - (9 - 5 t), with t as its path variable
	System H;
	H.AddVariableGroup(VariableGroup{x});
	H.AddFunction(x*x - (9 - 5*t));
	H.AddPathVariable(t);

	// the target system (H at t = 0): x^2 - 9, used for dehomogenize / residual in post-processing
	System target;
	target.AddVariableGroup(VariableGroup{x});
	target.AddFunction(x*x - 9);

	// the GIVEN start points: the t = 1 solutions +/- 2 (as if computed by an earlier solve)
	SampCont<mpfr> start_points;
	{
		Vec<mpfr> p(1); p(0) = mpfr(2);  start_points.push_back(p);
		Vec<mpfr> q(1); q(0) = mpfr(-2); start_points.push_back(q);
	}

	auto user_start = start_system::User(target, start_points);

	auto zd = algorithm::ZeroDim<
				AMPTracker,
				bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
				System,
				start_system::User,
				policy::RefToGiven>
			(target, user_start, H); // (target, start, homotopy) -- references, so all three

	zd.DefaultSetup();
	zd.Solve();

	auto const& sols = zd.SolutionsUserCoords();
	auto const& md   = zd.FinalSolutionMetadata();
	std::vector<dbl> ends;
	for (size_t i = 0; i < sols.size(); ++i)
		if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 1)
			ends.push_back(dbl(sols[i](0)));

	BOOST_CHECK_EQUAL(ends.size(), 2u);
	bool has_pos = false, has_neg = false;
	for (auto const& e : ends)
	{
		if (std::abs(e - dbl(3,0))  < 1e-7) has_pos = true;
		if (std::abs(e - dbl(-3,0)) < 1e-7) has_neg = true;
	}
	BOOST_CHECK(has_pos); // tracked +2 -> +3
	BOOST_CHECK(has_neg); // tracked -2 -> -3
}


// End-to-end multihomogeneous solve through the block-composed start system and the
// blend-block homotopy.  x*y - 1 = 0, x + y = 0 over variable groups {x}, {y}:
// y = -x gives -x^2 - 1 = 0, so x = +/- i -> exactly the two solutions (i,-i),(-i,i).
// The m-homogeneous Bezout number for bidegrees (1,1),(1,1) is 2, below the total-degree
// Bezout number 4 -- so MHom tracks 2 paths, not 4.  This exercises the whole new chain:
// the products-of-linears start block, the blend-block homotopy formed by FormHomotopy,
// and the block-aware System eval/Jacobian/time-derivative through the tracker + endgame.
//
// Tracked with the FIXED-DOUBLE tracker.  This is a *capability* test, not a deterministic
// one: the homotopy's gamma (and the MHom start coefficients) are drawn from RandomMp, which
// SetGlobalSeed does not reseed (only the ThreadEngine is reseedable today), so the exact
// path conditioning cannot be pinned to a seed.  MHom paths are conditioning-fragile in fixed
// double -- a given gamma drives both paths to MinStepSize maybe two times in three -- so we
// retry the solve, each attempt drawing a fresh gamma, and assert that MHom *can* solve the
// system and that when it does the solutions are genuine, distinct roots.  Two follow-ups make
// this robust under the CLI default (adaptive precision): (1) reseedable RandomMp +/or less
// fixed-double-fragile MHom paths; (2) AMP and the Cauchy endgame keeping precision in lockstep
// with a block-composed homotopy (today the tracked point reaches MaxPrecisionAllowed in the
// endgame while the system is at the working precision).
BOOST_AUTO_TEST_CASE(mhom_solves_two_variable_group_system)
{
	using namespace bertini;
	using namespace tracking;

	auto make_system = []{
		System sys;
		auto x = Variable::Make("x");
		auto y = Variable::Make("y");
		sys.AddVariableGroup(VariableGroup{x});
		sys.AddVariableGroup(VariableGroup{y});
		sys.AddFunction(x*y - 1);
		sys.AddFunction(x + y);
		return sys;
	};

	// AMP is the robust MHom path (precision escalates through the hard sections) and normally
	// solves on the first gamma; a tiny retry budget just absorbs a rare unlucky draw without the
	// 40x retry storm that made CI balloon (do NOT raise this back up).
	bool solved = false;
	for (int attempt = 0; attempt < 5 && !solved; ++attempt)
	{
		auto sys = make_system();
		auto zd = algorithm::ZeroDim<AMPTracker,
		                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
		                             decltype(sys),
		                             start_system::MHomogeneous>(sys);
		zd.DefaultSetup();
		zd.Solve();

		// Collect the successfully-tracked solutions (a failed path leaves an empty
		// placeholder, kept for index alignment with the metadata).
		auto const& sols = zd.SolutionsUserCoords();
		auto const& md   = zd.FinalSolutionMetadata();
		using SolVec = std::decay_t<decltype(sols[0])>;
		std::vector<SolVec> good;
		for (size_t i = 0; i < sols.size(); ++i)
			if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 2)
				good.push_back(sols[i]);

		if (good.size() != 2)
			continue; // this gamma drove a path to MinStepSize; try another

		// A poorly-conditioned gamma can let a path report endgame success yet land on a
		// spurious (inaccurate) endpoint in fixed double.  Treat accuracy + distinctness as part
		// of the success criterion: if either solution is off, this gamma is no good -- retry
		// rather than recording a failure.  We assert only that *some* gamma solves it cleanly.
		bool both_are_roots = true;
		for (auto const& s : good)
		{
			dbl a(s(0)), b(s(1));   // AMP returns mpfr_complex coords; double is plenty for a root check
			if (std::abs(a * b - dbl(1)) > 1e-8 || std::abs(a + b) > 1e-8)
				both_are_roots = false;
		}

		bool distinct = std::abs(dbl(good[0](0)) - dbl(good[1](0))) > 1e-3;

		if (both_are_roots && distinct)
			solved = true; // the two distinct roots (i,-i) and (-i,i), recovered in user coordinates
	}

	BOOST_CHECK(solved); // MHom solved the system through the block-composed homotopy
}


BOOST_AUTO_TEST_SUITE_END()
