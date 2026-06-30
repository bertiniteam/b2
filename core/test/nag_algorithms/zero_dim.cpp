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
#include "bertini2/trackers/observers.hpp"
#include "bertini2/trackers/events.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <typeindex>


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

	auto zd = algorithm::ZeroDimSolver<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);

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

	auto zd = algorithm::ZeroDimSolver<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, decltype(sys)>(sys);

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
	auto TD = start_system::TotalDegreeBinomial(sys);

	auto t = Variable::Make("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::HomotopySolver<
				TrackerT,
				bertini::endgame::EndgameSelector<TrackerT>::PSEG,
				decltype(sys)
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

	auto TD = start_system::TotalDegreeBinomial(sys);

	auto t = Variable::Make("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::HomotopySolver<
				TrackerT,
				bertini::endgame::EndgameSelector<TrackerT>::PSEG,
				decltype(sys)
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
// same engine, merely driven from start_system::User (start points come
// from the supplied list, not generated) and a user-supplied homotopy (taken as-is,
// not formed by homogenizing + coupling a start system).
//
// Parameter homotopy H(x,t) = x^2 - (9 - 5 t).  At t = 1 the roots are +/-2 (the given start
// points); at t = 0 they are +/-3 (the target x^2 - 9).  Tracking the two start points down to
// t = 0 must recover +/-3 -- i.e. ZeroDim moved the t=1 solutions to the t=0 parameter.
BOOST_AUTO_TEST_CASE(user_homotopy_parameter_homotopy_solves)
{
	using namespace bertini;
	using namespace tracking;
	using mpfr = bertini::complex_mp;

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

	auto zd = algorithm::HomotopySolver<
				AMPTracker,
				bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
				System>
			(target, user_start, H); // (target, start, homotopy) -- references, so all three

	zd.DefaultSetup();
	zd.Solve();

	auto const& sols = zd.SolutionsUserCoords();
	auto const& md   = zd.SolutionMetadata();
	std::vector<complex_dbl> ends;
	for (size_t i = 0; i < sols.size(); ++i)
		if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 1)
			ends.push_back(complex_dbl(sols[i](0)));

	BOOST_CHECK_EQUAL(ends.size(), 2u);
	bool has_pos = false, has_neg = false;
	for (auto const& e : ends)
	{
		if (std::abs(e - complex_dbl(3,0))  < 1e-7) has_pos = true;
		if (std::abs(e - complex_dbl(-3,0)) < 1e-7) has_neg = true;
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
// and the block-aware System eval/Jacobian/time-derivative through the AMP tracker + Cauchy endgame.
//
// Deterministic: SetGlobalSeed pins the homotopy gamma and the MHom start coefficients (RandomMp is
// reseedable now), so this is a single, reproducible solve -- no gamma-retry loop.
BOOST_AUTO_TEST_CASE(mhom_solves_two_variable_group_system)
{
	using namespace bertini;
	using namespace tracking;

	SetGlobalSeed(1); // reproducible homotopy gamma + MHom start coefficients

	System sys;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x + y);

	auto zd = algorithm::ZeroDimSolver<AMPTracker,
	                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
	                             decltype(sys)>(sys, bertini::start_system::MakeStartFactory<bertini::start_system::MHomogeneous>());
	zd.DefaultSetup();
	zd.Solve();

	// Collect the successfully-tracked solutions (a failed path leaves an empty placeholder).
	auto const& sols = zd.SolutionsUserCoords();
	auto const& md   = zd.SolutionMetadata();
	using SolVec = std::decay_t<decltype(sols[0])>;
	std::vector<SolVec> good;
	for (size_t i = 0; i < sols.size(); ++i)
		if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 2)
			good.push_back(sols[i]);

	BOOST_REQUIRE_EQUAL(good.size(), 2u); // both MHom paths solved (the m-homogeneous Bezout number)

	for (auto const& s : good)   // AMP returns complex_mp coords; double is plenty for a root check
	{
		complex_dbl a(s(0)), b(s(1));
		BOOST_CHECK_SMALL(std::abs(a * b - complex_dbl(1)), 1e-8);
		BOOST_CHECK_SMALL(std::abs(a + b), 1e-8);
	}
	BOOST_CHECK_GT(std::abs(complex_dbl(good[0](0)) - complex_dbl(good[1](0))), 1e-3); // the two distinct roots
}




// The decisive block-vs-function-tree check on the ACTUAL problematic homotopy: take the seed-6
// MHom blend homotopy (a BlendBlock of the homogenized target's polynomial block and the
// products-of-linears start), build its pure-function-tree twin via ExpandToFunctionTree(), and
// assert the two evaluate and (crucially) DIFFERENTIATE identically -- across random points and
// near t->0 (the endgame region), in double and mpfr.  If the block Jacobian were wrong (inflating
// ||J^{-1}|| and hence DigitsB), this is where it would show.  Together with the faithful
// ||J^{-1}|| estimate (amp_jacobian_estimate) and DegreeBound=2 (so Phi is tiny), agreement here
// means the DigitsB escalation reflects a GENUINE near-singular pass for that gamma, not a bug in
// the block representation.
BOOST_AUTO_TEST_CASE(mhom_homotopy_block_matches_function_tree)
{
	using namespace bertini;
	using namespace tracking;

	SetGlobalSeed(6); // the gamma that escalated to 140 digits in the probe

	System sys;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x + y);

	auto zd = algorithm::ZeroDimSolver<AMPTracker,
	                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
	                             decltype(sys)>(sys, bertini::start_system::MakeStartFactory<bertini::start_system::MHomogeneous>());
	zd.DefaultSetup();

	System const& H = zd.Homotopy();          // the block-composed blend homotopy
	System He = H.ExpandToFunctionTree();      // its pure-function-tree twin

	BOOST_CHECK_EQUAL(H.NumVariables(), He.NumVariables());
	BOOST_CHECK_EQUAL(H.NumTotalFunctions(), He.NumTotalFunctions());
	BOOST_CHECK_EQUAL(H.DegreeBound(), He.DegreeBound()); // drives Phi/Psi; must match

	auto compare = [](auto const& A, auto const& B, double tol)
	{
		BOOST_REQUIRE_EQUAL(A.rows(), B.rows());
		BOOST_REQUIRE_EQUAL(A.cols(), B.cols());
		for (Eigen::Index i = 0; i < A.rows(); ++i)
			for (Eigen::Index j = 0; j < A.cols(); ++j)
			{
				double err   = static_cast<double>(abs(A(i,j) - B(i,j)));
				double scale = 1.0 + static_cast<double>(abs(A(i,j)));
				BOOST_CHECK_SMALL(err / scale, tol);
			}
	};

	const int n = static_cast<int>(H.NumVariables());

	// double, including t very close to 0 (endgame region, where the spike lives)
	for (int trial = 0; trial < 10; ++trial)
	{
		Vec<complex_dbl> p = Vec<complex_dbl>::Random(n);
		complex_dbl t = (trial < 5) ? complex_dbl(0.4, -0.3) * complex_dbl(trial + 1)
		                    : complex_dbl(std::pow(10.0, -(trial - 1)), 0.0); // 1e-4 .. 1e-8
		compare(H.Eval(p, t),     He.Eval(p, t),     1e-11);
		compare(H.Jacobian(p, t), He.Jacobian(p, t), 1e-11);
	}

	// mpfr at 80 digits
	DefaultPrecision(80);
	H.precision(80);
	He.precision(80);
	for (int trial = 0; trial < 5; ++trial)
	{
		Vec<complex_mp> p = RandomOfUnits<complex_mp>(static_cast<unsigned int>(n));
		complex_mp t = RandomOfUnits<complex_mp>(1)(0);
		compare(H.Eval(p, t),     He.Eval(p, t),     1e-70);
		compare(H.Jacobian(p, t), He.Jacobian(p, t), 1e-70);
	}
}


// SolveReport: the end-of-solve diagnostic summary.  SummarizeSolve is pure aggregation, so unit-test
// it on synthetic metadata that includes a genuine tracking failure.
BOOST_AUTO_TEST_CASE(solve_report_buckets_metadata_and_flags_failures)
{
	using bertini::algorithm::SolutionMetaData;
	using bertini::algorithm::SummarizeSolve;
	using bertini::algorithm::MidpathCheckReport;
	using SC = bertini::SuccessCode;
	using CT = bertini::complex_dbl;

	auto set = [](SolutionMetaData<CT>& m, SC code, bool finite, int mult,
	              bool real, bool sing, double cond){
		m.endgame_success = code; m.is_finite = finite; m.multiplicity = mult;
		m.is_real = real; m.is_singular = sing; m.condition_number = cond; m.max_precision_used = 16;
	};

	std::vector<SolutionMetaData<CT>> md(8);
	set(md[0], SC::Success,            true,  1, true,  false, 1e3);   // finite, real
	set(md[1], SC::Success,            true,  1, false, false, 1e8);   // finite, complex
	set(md[2], SC::Success,            true,  2, false, true,  1e9);   // a double root: two paths...
	set(md[3], SC::Success,            true,  2, false, true,  1e9);   // ...counted once (sum 1/mult)
	set(md[4], SC::Success,            false, 1, false, false, 0);     // diverged (success + infinite)
	set(md[5], SC::GoingToInfinity,    false, 1, false, false, 0);     // diverged (clean)
	set(md[6], SC::MinStepSizeReached, false, 1, false, false, 0);     // FAILED -- a lost path
	set(md[7], SC::FailedToConverge,   false, 1, false, false, 0);     // FAILED

	auto r = SummarizeSolve(md, MidpathCheckReport{});
	BOOST_CHECK_EQUAL(r.num_paths_tracked,    8u);
	BOOST_CHECK_EQUAL(r.num_finite_endpoints, 4u);
	BOOST_CHECK_EQUAL(r.num_finite_solutions, 3u);   // 1 + 1 + (1/2 + 1/2)
	BOOST_CHECK_EQUAL(r.num_diverged,         2u);
	BOOST_CHECK_EQUAL(r.num_failed,           2u);
	BOOST_CHECK_EQUAL(r.num_real,             1u);
	BOOST_CHECK_EQUAL(r.num_singular,         2u);
	BOOST_CHECK_EQUAL(r.failures_by_reason[SC::MinStepSizeReached], 1u);
	BOOST_CHECK_EQUAL(r.failures_by_reason[SC::FailedToConverge],   1u);
	BOOST_CHECK(r.max_condition_number >= 1e9);
	BOOST_CHECK(!r.all_paths_resolved);              // failures present -> not trustworthy

	// drop the two failures: now every path resolved (and the midpath check passed by default)
	std::vector<SolutionMetaData<CT>> clean(md.begin(), md.begin() + 6);
	auto rc = SummarizeSolve(clean, MidpathCheckReport{});
	BOOST_CHECK_EQUAL(rc.num_failed, 0u);
	BOOST_CHECK(rc.all_paths_resolved);
}

// Regression: a path the security check truncated near infinity (SecurityMaxNormReached) is a
// divergence, not a failure.  With infinite-path truncation on by default, cyclic-n's infinite
// paths come back as SecurityMaxNormReached; before this fix they were bucketed as num_failed,
// which (wrongly) cleared all_paths_resolved and broke the cyclic-5 example's `num_failed == 0` gate.
BOOST_AUTO_TEST_CASE(solve_report_security_truncation_counts_as_diverged)
{
	using bertini::algorithm::SolutionMetaData;
	using bertini::algorithm::SummarizeSolve;
	using bertini::algorithm::MidpathCheckReport;
	using SC = bertini::SuccessCode;
	using CT = bertini::complex_dbl;

	auto set = [](SolutionMetaData<CT>& m, SC code, bool finite){
		m.endgame_success = code; m.is_finite = finite; m.multiplicity = 1;
		m.is_real = false; m.is_singular = false; m.condition_number = 1e3; m.max_precision_used = 16;
	};

	std::vector<SolutionMetaData<CT>> md(4);
	set(md[0], SC::Success,                true);   // finite
	set(md[1], SC::GoingToInfinity,        false);  // diverged (clean)
	set(md[2], SC::SecurityMaxNormReached, false);  // diverged (security-truncated near infinity)
	set(md[3], SC::SecurityMaxNormReached, false);  // diverged (security-truncated near infinity)

	auto r = SummarizeSolve(md, MidpathCheckReport{});
	BOOST_CHECK_EQUAL(r.num_finite_endpoints, 1u);
	BOOST_CHECK_EQUAL(r.num_diverged,         3u);   // 1 GoingToInfinity + 2 SecurityMaxNormReached
	BOOST_CHECK_EQUAL(r.num_failed,           0u);   // truncations are NOT failures
	BOOST_CHECK(r.failures_by_reason.find(SC::SecurityMaxNormReached) == r.failures_by_reason.end());
	BOOST_CHECK(r.all_paths_resolved);
}

// A real solve: x^2-1, y^2-1 -> exactly 4 finite roots, no losses.
BOOST_AUTO_TEST_CASE(solve_report_from_a_real_solve)
{
	using namespace bertini;

	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddFunction(pow(y, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto const& meta = zd.SolutionMetadata();
	auto r = zd.Report();

	BOOST_CHECK_EQUAL(r.num_paths_tracked, meta.size());                       // total degree 2*2 = 4
	BOOST_CHECK_EQUAL(r.num_finite_endpoints + r.num_diverged + r.num_failed,  // the buckets partition
	                  r.num_paths_tracked);
	BOOST_CHECK_EQUAL(r.num_finite_solutions, 4u);
	BOOST_CHECK_EQUAL(r.num_failed, 0u);
	BOOST_CHECK(r.all_paths_resolved);

	std::ostringstream oss; oss << r;                                          // it prints
	BOOST_CHECK(!oss.str().empty());
}

// ADR-0038 regression: AMP Criterion B in the tracker cost model must use the latest Newton residual
// (norm_delta_z), NOT size_proportion.  size_proportion = err_est/|dt|^(p+1) blows up to ~1e51 in the
// endgame roundoff regime; passed raw into Criterion B it added ~25 spurious digits and forced ~39 of
// cyclic5's 70 finite paths into multiprecision even though all 70 converge in pure double.  After the
// fix only a couple paths need MP.  Guard: solve cyclic5 with AMP and assert almost every path stays
// in double.
BOOST_AUTO_TEST_CASE(cyclic5_amp_does_not_overescalate_precision)
{
	using namespace bertini;
	const int N = 5;

	VariableGroup x;
	for (int i = 0; i < N; ++i)
		x.push_back(node::Variable::Make("x" + std::to_string(i)));

	System sys;
	// f_k = sum_i prod_{j=0}^{k-1} x_{(i+j) mod N},  for k = 1 .. N-1
	for (int k = 1; k < N; ++k)
	{
		std::shared_ptr<node::Node> f;
		for (int i = 0; i < N; ++i)
		{
			std::shared_ptr<node::Node> term = x[i];
			for (int j = 1; j < k; ++j)
				term = term * x[(i + j) % N];
			f = (i == 0) ? term : (f + term);
		}
		sys.AddFunction(f);
	}
	// closing relation: (prod_i x_i) - 1
	std::shared_ptr<node::Node> prod = x[0];
	for (int i = 1; i < N; ++i)
		prod = prod * x[i];
	sys.AddFunction(prod - 1);
	sys.AddVariableGroup(x);

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	zd.Solve();

	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size(), 70u);

	unsigned escalated = 0;
	for (auto const& m : zd.SolutionMetadata())
		if (m.precision_changed)
			++escalated;
	BOOST_TEST_MESSAGE("cyclic5 AMP finite paths that escalated to multiprecision: " << escalated << " / 70");
	BOOST_CHECK_LT(escalated, 10u);   // ~39 before the ADR-0038 fix; ~2 after
}

// The filtered-solution convenience accessors, on the same x^2-1, y^2-1 solve: all four roots
// (+-1, +-1) are finite, real, and nonsingular.
BOOST_AUTO_TEST_CASE(filtered_solution_accessors)
{
	using namespace bertini;

	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddFunction(pow(y, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto r = zd.Report();
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size(),      4u);
	BOOST_CHECK_EQUAL(zd.RealSolutions().size(),        4u);
	BOOST_CHECK_EQUAL(zd.NonsingularSolutions().size(), 4u);
	BOOST_CHECK_EQUAL(zd.SingularSolutions().size(),    0u);

	// the accessors agree with the report's counts
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size(), r.num_finite_endpoints);
	BOOST_CHECK_EQUAL(zd.RealSolutions().size(),   r.num_real);
	BOOST_CHECK_EQUAL(zd.SingularSolutions().size(), r.num_singular);

	// singular + nonsingular partition the finite set
	BOOST_CHECK_EQUAL(zd.SingularSolutions().size() + zd.NonsingularSolutions().size(),
	                  zd.FiniteSolutions().size());

	// user vs internal coordinates: same count, different representation
	BOOST_CHECK_EQUAL(zd.FiniteSolutions(false).size(), zd.FiniteSolutions(true).size());
}

// The Bertini 1.7-compatible solution-file writers (output::Classic): each file is count-led and
// machine-readable -- first line is the solution count, then one block of NumVariables "re im"
// coordinate lines per solution.  On x^2-1, y^2-1 the four roots (+-1,+-1) are all finite, real,
// and nonsingular.
BOOST_AUTO_TEST_CASE(classic_solution_file_output)
{
	using namespace bertini;
	namespace out = bertini::algorithm::output;

	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddFunction(pow(y, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	zd.Solve();

	// first whitespace-delimited token is the solution count
	auto first_count = [](std::string const& s){ std::istringstream iss(s); long n=-1; iss >> n; return n; };
	// "re im" coordinate lines contain a space; count lines, path-index lines, and blanks do not
	auto coord_lines = [](std::string const& s){
		std::istringstream iss(s); std::string line; std::size_t c=0;
		while (std::getline(iss, line)) if (line.find(' ') != std::string::npos) ++c;
		return c; };

	std::ostringstream fin, real, nonsing, sing, raw;
	out::Classic<decltype(zd)>::FiniteSolutions(fin, zd);
	out::Classic<decltype(zd)>::RealFiniteSolutions(real, zd);
	out::Classic<decltype(zd)>::NonsingularSolutions(nonsing, zd);
	out::Classic<decltype(zd)>::SingularSolutions(sing, zd);
	out::Classic<decltype(zd)>::RawSolutions(raw, zd);

	BOOST_CHECK_EQUAL(first_count(fin.str()),     4);
	BOOST_CHECK_EQUAL(first_count(real.str()),    4);
	BOOST_CHECK_EQUAL(first_count(nonsing.str()), 4);
	BOOST_CHECK_EQUAL(first_count(sing.str()),    0);
	BOOST_CHECK_EQUAL(first_count(raw.str()),     4);

	// the count in the file agrees with the corresponding accessor
	BOOST_CHECK_EQUAL(first_count(fin.str()),  static_cast<long>(zd.FiniteSolutions().size()));
	BOOST_CHECK_EQUAL(first_count(sing.str()), static_cast<long>(zd.SingularSolutions().size()));

	// exactly one coordinate block (NumVariables lines) per finite solution
	BOOST_CHECK_EQUAL(coord_lines(fin.str()), sys.NumVariables() * 4u);
}

// InfiniteSolutions: the at-infinity complement of FiniteSolutions.  The system {x*y - 1, x - 1}
// has Bezout number 2 but exactly one affine solution (1,1); the second total-degree path must
// diverge -- so there is one finite endpoint and one at infinity.
BOOST_AUTO_TEST_CASE(infinite_solutions_at_infinity)
{
	using namespace bertini;

	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x - 1);
	sys.AddVariableGroup(VariableGroup{x, y});

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	// SecurityLevel <= 0 truncates paths heading to infinity as failures (SecurityMaxNormReached);
	// this test wants the at-infinity path actually computed, so raise the level to keep tracking it
	// to its (infinite) endpoint, where it is classified as a divergence.
	endgame::SecurityConfig sec;
	sec.level = 1;
	zd.GetEndgame().Set(sec);
	zd.Solve();

	auto r = zd.Report();
	BOOST_CHECK_EQUAL(zd.SolutionMetadata().size(), 2u);   // Bezout 2
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size(),       1u);
	BOOST_CHECK_EQUAL(zd.InfiniteSolutions().size(),     1u);

	// finite + infinite partition the whole list (no failed paths on this clean solve)
	BOOST_CHECK_EQUAL(r.num_failed, 0u);
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size() + zd.InfiniteSolutions().size(),
	                  zd.SolutionMetadata().size());

	// the at-infinity count matches the report's diverged bucket
	BOOST_CHECK_EQUAL(zd.InfiniteSolutions().size(), r.num_diverged);

	// user vs internal coordinates: same count
	BOOST_CHECK_EQUAL(zd.InfiniteSolutions(false).size(), zd.InfiniteSolutions(true).size());
}

// multiplicity_representative: a multiplicity-m solution arrives as m coincident endpoints, and
// the solver must mark exactly ONE of them the representative (the rest false), so a consumer can
// collapse the cluster to one row.  {x^2, y^2} has a single solution (0,0) of multiplicity 4.
BOOST_AUTO_TEST_CASE(multiplicity_representative_marks_one_per_cluster)
{
	using namespace bertini;

	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2));
	sys.AddFunction(pow(y, 2));
	sys.AddVariableGroup(VariableGroup{x, y});

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto const& md = zd.SolutionMetadata();
	BOOST_CHECK_EQUAL(md.size(), 4u);                                // Bezout 2*2 = 4 paths

	unsigned representatives = 0, duplicates = 0;
	for (auto const& m : md)
	{
		if (!m.is_finite) continue;
		if (m.multiplicity_representative) { ++representatives; BOOST_CHECK_EQUAL(m.multiplicity, 4); }
		else                                 ++duplicates;
	}
	BOOST_CHECK_EQUAL(representatives, 1u);                          // exactly one representative
	BOOST_CHECK_EQUAL(duplicates,     3u);                          // the other m-1 copies

	// the representative count equals the number of DISTINCT finite solutions in the report
	BOOST_CHECK_EQUAL(representatives, zd.Report().num_finite_solutions);

	// a clean simple-root solve marks every finite endpoint a representative (nothing to merge):
	// x^2 - 1, y^2 - 1 has four distinct simple roots.
	System simple;
	simple.AddFunction(pow(x, 2) - 1);
	simple.AddFunction(pow(y, 2) - 1);
	simple.AddVariableGroup(VariableGroup{x, y});
	auto zd2 = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(simple);
	zd2.DefaultSetup();
	zd2.Solve();
	unsigned reps2 = 0;
	for (auto const& m : zd2.SolutionMetadata())
		if (m.is_finite && m.multiplicity_representative) ++reps2;
	BOOST_CHECK_EQUAL(reps2, 4u);
}

// Regression: a fixed-multiple (MultiplePrecisionTracker) zero-dim solve used to throw at the start
// of tracking -- "start point ... differing precision from default (20!=16)" -- because the tracker
// (built at DefaultPrecision) and the config-driven ambient/thread precision (DoublePrecision)
// disagreed.  ZeroDimConfig::initial_ambient_precision now defaults to the multiprecision
// DefaultPrecision, so the precision is uniform everywhere and the solve completes -- at whatever
// precision is the default when the solver is constructed.
BOOST_AUTO_TEST_CASE(fixed_multiple_precision_solves_uniformly)
{
	using namespace bertini;
	using MPTracker = tracking::MultiplePrecisionTracker;

	auto solve_at = [](unsigned at_precision) -> unsigned long long {
		auto saved = DefaultPrecision();
		DefaultPrecision(at_precision);
		auto x = node::Variable::Make("x");
		System sys;
		sys.AddFunction(pow(x, 2) - 1);
		sys.AddVariableGroup(VariableGroup{x});
		auto zd = algorithm::ZeroDimSolver<MPTracker, endgame::EndgameSelector<MPTracker>::Cauchy,
		                             System>(sys);
		zd.DefaultSetup();
		zd.Solve();                                  // used to throw on the precision mismatch
		auto n = zd.Report().num_finite_solutions;
		DefaultPrecision(saved);
		return n;
	};

	BOOST_CHECK_EQUAL(solve_at(30), 2u);             // a typical multiprecision default
	BOOST_CHECK_EQUAL(solve_at(60), 2u);             // and a higher one -- still uniform, still solves
}

// FixedPrecisionConfig.precision sets a fixed-multiple solve's precision *after* construction -- no
// DefaultPrecision()-before-construct dance.  This is the "real precision setter" ADR-0030 deferred:
// the algorithm lifts the tracker, the systems, and the ambient/start-point precision to the config's
// value at PreSolveSetup.
BOOST_AUTO_TEST_CASE(fixed_multiple_precision_set_via_config)
{
	using namespace bertini;
	using MPTracker = tracking::MultiplePrecisionTracker;

	auto saved = DefaultPrecision();
	DefaultPrecision(30);                            // construct at a modest default
	auto x = node::Variable::Make("x");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x});
	auto zd = algorithm::ZeroDimSolver<MPTracker, endgame::EndgameSelector<MPTracker>::Cauchy,
	                             System>(sys);
	zd.DefaultSetup();

	// the tracker reports its precision honestly
	BOOST_CHECK_EQUAL(zd.GetTracker().template Get<tracking::FixedPrecisionConfig>().precision, 30u);

	// choose a different fixed precision via the config, then solve at it
	auto fp = zd.GetTracker().template Get<tracking::FixedPrecisionConfig>();
	fp.precision = 70;
	zd.GetTracker().template Set<tracking::FixedPrecisionConfig>(fp);

	zd.Solve();
	BOOST_CHECK_EQUAL(zd.Report().num_finite_solutions, 2u);
	BOOST_CHECK_EQUAL(zd.GetTracker().CurrentPrecision(), 70u);   // the solve ran at 70
	DefaultPrecision(saved);
}

// Double-precision tracking is fixed at DoublePrecision(); asking the config for any other precision
// is rejected (use mptype multiple/adaptive for more digits).
BOOST_AUTO_TEST_CASE(double_precision_rejects_a_different_precision)
{
	using namespace bertini;
	auto x = node::Variable::Make("x");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x});
	tracking::DoublePrecisionTracker tracker(sys);

	BOOST_CHECK_EQUAL(tracker.template Get<tracking::FixedPrecisionConfig>().precision, DoublePrecision());
	tracking::FixedPrecisionConfig fp;
	fp.precision = 50;
	BOOST_CHECK_THROW(tracker.PrecisionSetup(fp), std::runtime_error);
	fp.precision = DoublePrecision();
	BOOST_CHECK_NO_THROW(tracker.PrecisionSetup(fp));   // the native precision is fine
}

// PROBE observer (branch perf/amp-block-precision-escalation): record, per successful step, the
// |t|, working precision, and the tracker's condition-number estimate -- so we can SEE whether the
// condition number (||J|| * ||J^{-1}||) spikes then RECOVERS along the actual seed-6 path, and at
// what |t| (mid-path vs the t->0 endgame region).
template <class TrackerT>






BOOST_AUTO_TEST_SUITE_END()



// Path-crossing (midpath) detection and the re-track machinery (RunMidpathResolution).
//
// Two distinct paths landing on the same point at the endgame boundary is a probability-0 event:
// it signals a path crossing (under-resolved tracking), which the algorithm is supposed to detect
// and re-track.  These tests guard the two historical defects in that machinery:
//   * Bug 1: MidpathChecker computed `same_start` with an inverted comparison, so a genuine
//     crossing (distinct starts, coincident boundary points) was flagged rerun==false and never
//     re-tracked.
//   * Bug 2: MidpathChecker::Check did not reset its pass flag / crossed-path list between calls,
//     so the resolve loop never saw a clean pass and operated on stale data.
// The first three cases are deterministic (hand-built boundary + start data) so they always pass
// once the logic is correct, independent of tracker numerics.  The end-to-end "provoke a spurious
// crossing and watch it get resolved" scenario is exercised in the Python suite on the cyclic
// system (see python/test/zero_dim/crossed_paths_test.py), where a cyclic builder exists and the
// phenomenon is documented (ADR-0017); the C++ wiring case below just confirms a clean solve
// reports no crossings and the new accessors work.
BOOST_AUTO_TEST_SUITE(crossed_paths)

using bertini::complex_dbl;
using bertini::Vec;
using bertini::SuccessCode;
using MidPathConfig = bertini::algorithm::MidPathConfig;
using BoundaryMD = bertini::algorithm::EGBoundaryMetaData<complex_dbl>;
using Checker = bertini::algorithm::MidpathChecker<double, complex_dbl, BoundaryMD>;

namespace {
	// Minimal stand-in providing just the StartPoint<ComplexT>(index) interface MidpathChecker
	// needs, so the test controls both the boundary points and the start points.
	struct MockStartSystem
	{
		std::vector<Vec<complex_dbl>> pts;

		template<typename ComplexT>
		Vec<ComplexT> StartPoint(unsigned long long i) const
		{
			return pts[i].template cast<ComplexT>();
		}
	};

	Vec<complex_dbl> Pt(complex_dbl a, complex_dbl b)
	{
		Vec<complex_dbl> v(2);
		v << a, b;
		return v;
	}

	BoundaryMD MakeBoundaryPoint(Vec<complex_dbl> const& p)
	{
		return BoundaryMD(p, SuccessCode::Success, 0.01, 16);
	}
}


// Bug 1: two coincident boundary points from DISTINCT start points is a genuine crossing -- it must
// be detected (Check fails) and both involved paths must be flagged for re-tracking (rerun==true).
BOOST_AUTO_TEST_CASE(distinct_starts_coincident_endpoints_are_rerun)
{
	Checker mp{MidPathConfig()};

	// paths 0 and 1 land on (essentially) the same boundary point; path 2 is elsewhere.
	auto p = Pt(1.0, 1.0);
	auto p_nudged = Pt(1.0 + 1e-9, 1.0 - 1e-9);
	auto p_far = Pt(5.0, 5.0);
	Checker::BoundaryData boundary{MakeBoundaryPoint(p), MakeBoundaryPoint(p_nudged), MakeBoundaryPoint(p_far)};

	MockStartSystem starts{{Pt(1.0, 0.0), Pt(2.0, 0.0), Pt(3.0, 0.0)}}; // all distinct

	bool passed = mp.Check(boundary, starts);

	BOOST_CHECK(!passed);
	auto crossed = mp.GetCrossedPaths();
	BOOST_REQUIRE_EQUAL(crossed.size(), 2u); // paths 0 and 1
	for (auto const& c : crossed)
	{
		BOOST_CHECK(c.index() == 0 || c.index() == 1);
		BOOST_CHECK(c.rerun()); // distinct starts -> genuine crossing -> must re-track
	}
}


// Two coincident boundary points that began at the SAME start point are not a crossing to re-track:
// the crossing is still detected (Check fails) but rerun must be false.
BOOST_AUTO_TEST_CASE(same_start_coincident_endpoints_are_not_rerun)
{
	Checker mp{MidPathConfig()};

	auto p = Pt(1.0, 1.0);
	auto p_nudged = Pt(1.0 + 1e-9, 1.0 - 1e-9);
	auto p_far = Pt(5.0, 5.0);
	Checker::BoundaryData boundary{MakeBoundaryPoint(p), MakeBoundaryPoint(p_nudged), MakeBoundaryPoint(p_far)};

	// paths 0 and 1 share a start point.
	MockStartSystem starts{{Pt(1.0, 0.0), Pt(1.0, 0.0), Pt(3.0, 0.0)}};

	bool passed = mp.Check(boundary, starts);

	BOOST_CHECK(!passed);
	auto crossed = mp.GetCrossedPaths();
	BOOST_REQUIRE_EQUAL(crossed.size(), 2u);
	for (auto const& c : crossed)
		BOOST_CHECK(!c.rerun()); // same start -> not a re-trackable crossing
}


// Bug 2: Check must reset its state each call.  A first call that finds a crossing must not leave
// the checker permanently "failed": a subsequent call on clean data must pass and report no
// crossings.
BOOST_AUTO_TEST_CASE(check_resets_state_between_calls)
{
	Checker mp{MidPathConfig()};

	auto p = Pt(1.0, 1.0);
	auto p_nudged = Pt(1.0 + 1e-9, 1.0 - 1e-9);
	auto p_far = Pt(5.0, 5.0);
	MockStartSystem starts{{Pt(1.0, 0.0), Pt(2.0, 0.0), Pt(3.0, 0.0)}};

	Checker::BoundaryData crossing{MakeBoundaryPoint(p), MakeBoundaryPoint(p_nudged), MakeBoundaryPoint(p_far)};
	BOOST_CHECK(!mp.Check(crossing, starts));
	BOOST_CHECK(!mp.GetCrossedPaths().empty());

	// now feed clean (all-distinct) data: must pass, and the stale crossing list must be gone.
	Checker::BoundaryData clean{MakeBoundaryPoint(Pt(1.0, 1.0)), MakeBoundaryPoint(Pt(5.0, 5.0)), MakeBoundaryPoint(Pt(9.0, 9.0))};
	BOOST_CHECK(mp.Check(clean, starts));
	BOOST_CHECK(mp.GetCrossedPaths().empty());
}


// Wiring: a clean solve populates the midpath report (no crossings) and the boundary accessors work.
BOOST_AUTO_TEST_CASE(clean_solve_reports_no_crossings)
{
	using namespace bertini;
	using TrackerT = tracking::DoublePrecisionTracker;

	auto sys = system::Precon::GriewankOsborn();
	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto const& report = zd.EndgameBoundaryMetadata();
	BOOST_CHECK(report.passed);
	BOOST_CHECK_EQUAL(report.num_crossings_detected, 0u);
	BOOST_CHECK_EQUAL(report.num_resolve_attempts, 0u);
	BOOST_CHECK(zd.EndgameBoundarySolutions().size() > 0u);
}

BOOST_AUTO_TEST_SUITE_END()



// Solution-metadata classification (is_finite / is_real / is_singular / multiplicity) and the
// PostProcessing config knobs that drive it.  These pin the Bertini-1 behaviour: an endpoint is
// at infinity if the infinity norm of its dehomogenized coordinates exceeds
// endpoint_finite_threshold; real if the imaginary parts are below real_threshold; singular if it
// is a multiple endpoint or its condition number exceeds condition_number_threshold; and two
// endpoints are the same when their dehomogenized coordinates agree to
// final_tolerance * same_point_tolerance_multiplier (infinity norm).
BOOST_AUTO_TEST_SUITE(zero_dim_solution_metadata)

using TrackerT = bertini::tracking::DoublePrecisionTracker;
using PostProcessing = bertini::algorithm::PostProcessingConfig;

// tally the classification flags over the successfully-tracked endpoints
struct Counts { int success = 0, finite = 0, real = 0, singular = 0; };
template<typename MDVec>
Counts Tally(MDVec const& md)
{
	Counts c;
	for (auto const& m : md)
	{
		if (m.endgame_success != bertini::SuccessCode::Success) continue;
		++c.success;
		if (m.is_finite)   ++c.finite;
		if (m.is_real)     ++c.real;
		if (m.is_singular) ++c.singular;
	}
	return c;
}

// x^2 - 1 -> roots +/-1: two finite, real, nonsingular solutions.
BOOST_AUTO_TEST_CASE(finite_real_nonsingular)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x - 1);

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto c = Tally(zd.SolutionMetadata());
	BOOST_CHECK_EQUAL(c.success, 2);
	BOOST_CHECK_EQUAL(c.finite, 2);
	BOOST_CHECK_EQUAL(c.real, 2);
	BOOST_CHECK_EQUAL(c.singular, 0);
	for (auto const& m : zd.SolutionMetadata())
		if (m.endgame_success == SuccessCode::Success)
			BOOST_CHECK_EQUAL(m.multiplicity, 1);
}


// x^2 + 1 -> roots +/-i: two finite, NON-real, nonsingular solutions.
BOOST_AUTO_TEST_CASE(finite_complex_not_real)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x + 1);

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto c = Tally(zd.SolutionMetadata());
	BOOST_CHECK_EQUAL(c.success, 2);
	BOOST_CHECK_EQUAL(c.finite, 2);
	BOOST_CHECK_EQUAL(c.real, 0);     // +/- i are not real
	BOOST_CHECK_EQUAL(c.singular, 0);
}


// x^2 = 0 -> a double root at 0: singular (multiplicity 2, and ill-conditioned), finite, real.
BOOST_AUTO_TEST_CASE(singular_double_root)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x);

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto const& md = zd.SolutionMetadata();
	auto c = Tally(md);
	BOOST_CHECK(c.success >= 1);                 // the endgame should reach the singular endpoint
	BOOST_CHECK_EQUAL(c.singular, c.success);    // every successful endpoint here is singular
	BOOST_CHECK_EQUAL(c.finite, c.success);      // ... and finite (at 0)
	if (c.success == 2)                          // both paths clustered -> multiplicity 2
		for (auto const& m : md)
			if (m.endgame_success == SuccessCode::Success)
				BOOST_CHECK_EQUAL(m.multiplicity, 2);
}


// endpoint_finite_threshold is actually applied: lower it below the solutions' norm and the
// finite roots get classified as at infinity.
BOOST_AUTO_TEST_CASE(endpoint_finite_threshold_is_applied)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x - 1);                    // roots +/-1, infinity norm 1

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	auto pp = zd.Get<PostProcessing>();
	pp.endpoint_finite_threshold = 0.5;          // 1 > 0.5 -> "at infinity"
	zd.Set(pp);
	zd.Solve();

	auto c = Tally(zd.SolutionMetadata());
	BOOST_CHECK_EQUAL(c.success, 2);
	BOOST_CHECK_EQUAL(c.finite, 0);              // the lowered threshold reclassifies both as infinite
}


// condition_number_threshold is actually applied: lower it below the (well-conditioned) roots'
// condition number and they get classified as singular.
BOOST_AUTO_TEST_CASE(condition_number_threshold_is_applied)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x - 1);                    // simple, well-conditioned roots

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();
	auto pp = zd.Get<PostProcessing>();
	pp.condition_number_threshold = 1e-3;        // any condition number exceeds this
	zd.Set(pp);
	zd.Solve();

	auto c = Tally(zd.SolutionMetadata());
	BOOST_CHECK_EQUAL(c.success, 2);
	BOOST_CHECK_EQUAL(c.singular, 2);            // reclassified singular purely by the lowered threshold
}


// the PostProcessing config round-trips through Get/Set.
BOOST_AUTO_TEST_CASE(postprocessing_config_roundtrip)
{
	using namespace bertini;
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(x*x - 1);

	auto zd = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys)>(sys);
	zd.DefaultSetup();

	auto pp = zd.Get<PostProcessing>();
	pp.endpoint_finite_threshold       = 12345.0;
	pp.same_point_tolerance_multiplier = 7.0;
	pp.condition_number_threshold      = 99.0;
	pp.real_threshold                  = 1e-3;
	zd.Set(pp);

	auto pp2 = zd.Get<PostProcessing>();
	BOOST_CHECK_CLOSE(pp2.endpoint_finite_threshold,       12345.0, 1e-10);
	BOOST_CHECK_CLOSE(pp2.same_point_tolerance_multiplier, 7.0,     1e-10);
	BOOST_CHECK_CLOSE(pp2.condition_number_threshold,      99.0,    1e-10);
	BOOST_CHECK_CLOSE(pp2.real_threshold,                  1e-3,    1e-10);
}


// the PostProcessing defaults match Bertini 1 (guards against the inverted endpoint-finite default
// that used to ship, 1e-5 instead of 1e5).
BOOST_AUTO_TEST_CASE(postprocessing_config_defaults_match_bertini1)
{
	bertini::algorithm::PostProcessingConfig pp;
	BOOST_CHECK_CLOSE(pp.endpoint_finite_threshold,       1e5,  1e-8);
	BOOST_CHECK_CLOSE(pp.same_point_tolerance_multiplier, 10.0, 1e-8);
	BOOST_CHECK_CLOSE(pp.condition_number_threshold,      1e8,  1e-8);
	BOOST_CHECK_CLOSE(pp.real_threshold,                  1e-8, 1e-8);
}


// Observe a whole zero-dim solve: AlgorithmStarted once, a PathStarted/PathComplete
// pair per path, AlgorithmComplete once.  The observer attaches to the ZeroDim itself
// (the AnyZeroDim emitter), not to the tracker.
namespace {

struct ZeroDimLifecycleCounter : public bertini::Observer<bertini::algorithm::AnyZeroDim>
{
	int started = 0, completed = 0, path_begin = 0, path_end = 0;

	bertini::ObserveResult Observe(bertini::AnyEvent const& e) override
	{
		using namespace bertini::algorithm;
		if (dynamic_cast<const AlgorithmStarted<AnyZeroDim>*>(&e))        ++started;
		else if (dynamic_cast<const AlgorithmComplete<AnyZeroDim>*>(&e))  ++completed;
		else if (dynamic_cast<const PathStarted<AnyZeroDim>*>(&e))      ++path_begin;
		else if (dynamic_cast<const PathComplete<AnyZeroDim>*>(&e))       ++path_end;
		return bertini::ObserveResult::KeepObserving;
	}
};

} // anon namespace

BOOST_AUTO_TEST_CASE(zerodim_emits_lifecycle_events)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	auto zd = algorithm::ZeroDimSolver<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy,
	                             decltype(sys)>(sys);
	zd.DefaultSetup();

	ZeroDimLifecycleCounter counter;
	zd.AddObserver(counter);   // attaches to the ZeroDim (accepts an AnyZeroDim observer)

	zd.Solve();

	BOOST_CHECK_EQUAL(counter.started,   1);
	BOOST_CHECK_EQUAL(counter.completed, 1);
	BOOST_CHECK_GT(counter.path_begin, 0);
	BOOST_CHECK_EQUAL(counter.path_begin, counter.path_end);  // every path that began also completed
}


BOOST_AUTO_TEST_CASE(zerodim_rejects_a_tracker_observer)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	auto zd = algorithm::ZeroDimSolver<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy,
	                             decltype(sys)>(sys);

	// a tracker observer is for a tracker, not the ZeroDim -> rejected
	GoryDetailLogger<TrackerT> tracker_observer;
	BOOST_CHECK_THROW(zd.AddObserver(tracker_observer), bertini::IncompatibleObserver);
}


// --- ZeroDimSolver feasibility + square-up behaviors (the new in-scope features) ----------------

// A pleasant SQUARE system solves directly, with no randomization.
BOOST_AUTO_TEST_CASE(square_system_solves_without_randomization)
{
	using namespace bertini;
	using namespace tracking;

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y - 1);   // unit circle
	sys.AddFunction(x - y);           // meets the line x = y in two points

	auto zd = algorithm::ZeroDimSolver<AMPTracker, endgame::EndgameSelector<AMPTracker>::Cauchy, System>(sys);
	BOOST_CHECK(!zd.WasRandomized());
	zd.Solve();
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size(), 2u);
}

// An OVER-determined system is squared up by randomization, and the extraneous solutions the
// squaring introduces are filtered out -- finite_solutions returns only the genuine roots.
BOOST_AUTO_TEST_CASE(overdetermined_system_is_squared_and_filtered)
{
	using namespace bertini;
	using namespace tracking;

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y - 2);   // three equations, two variables -> over-determined
	sys.AddFunction(x - y);
	sys.AddFunction(x*y - 1);
	// the only common (genuine) solutions are (1,1) and (-1,-1).

	auto zd = algorithm::ZeroDimSolver<AMPTracker, endgame::EndgameSelector<AMPTracker>::Cauchy, System>(sys);
	BOOST_CHECK(zd.WasRandomized());
	zd.Solve();

	// the randomized (square) system had MORE endpoints than genuine solutions ...
	BOOST_CHECK_GT(zd.SolutionsUserCoords().size(), 2u);
	// ... but after filtering against the original system, exactly the two genuine roots remain.
	auto fin = zd.FiniteSolutions();
	BOOST_CHECK_EQUAL(fin.size(), 2u);
	bool has_pp = false, has_nn = false;
	for (auto const& p : fin)
	{
		auto a = complex_dbl(p(0)), b = complex_dbl(p(1));
		if (std::abs(a - complex_dbl(1))  < 1e-6 && std::abs(b - complex_dbl(1))  < 1e-6) has_pp = true;
		if (std::abs(a - complex_dbl(-1)) < 1e-6 && std::abs(b - complex_dbl(-1)) < 1e-6) has_nn = true;
	}
	BOOST_CHECK(has_pp);
	BOOST_CHECK(has_nn);

	// the filtered-out points are flagged is_nonsolution (not is_finite=false) and surfaced by
	// Nonsolutions(); they stay geometrically finite.  How many of the squaring's extra roots land
	// finite vs. diverge is RNG-dependent, so assert the ROBUST invariants, not an exact count:
	unsigned num_nonsol = 0;
	for (auto const& m : zd.SolutionMetadata())
		if (m.is_nonsolution) { ++num_nonsol; BOOST_CHECK(m.is_finite); }   // a nonsolution is finite
	BOOST_CHECK_EQUAL(num_nonsol, zd.Nonsolutions().size());
	BOOST_CHECK_EQUAL(zd.Report().num_nonsolutions, num_nonsol);
	// every endpoint is exactly one of: a genuine finite solution, a finite nonsolution, or at
	// infinity -- so the three accessors partition the full per-path list.
	BOOST_CHECK_EQUAL(zd.FiniteSolutions().size() + zd.Nonsolutions().size() + zd.InfiniteSolutions().size(),
	                  zd.SolutionsUserCoords().size());
}

// An UNDER-determined system has a positive-dimensional solution set, so ZeroDimSolver refuses it
// at construction with a helpful error (ConsistencyCheck).
BOOST_AUTO_TEST_CASE(underdetermined_system_raises)
{
	using namespace bertini;
	using namespace tracking;

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y - 1);   // one equation, two variables -> under-determined

	BOOST_CHECK_THROW(
		(algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys)),
		std::runtime_error);
}

// A system that is SQUARE by equation count but whose variety is still positive-dimensional (here a
// repeated equation) is rejected by the generic-point Jacobian rank check.
BOOST_AUTO_TEST_CASE(positive_dimensional_square_system_raises)
{
	using namespace bertini;
	using namespace tracking;

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*y);   // two copies: square count (2 eqns, 2 vars) but {xy = 0} is a curve
	sys.AddFunction(x*y);

	BOOST_CHECK_THROW(
		(algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>(sys)),
		std::runtime_error);
}


BOOST_AUTO_TEST_SUITE_END()
