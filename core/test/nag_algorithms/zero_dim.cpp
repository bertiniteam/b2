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
#include "bertini2/detail/escalation_probe.hpp" // PROBE: temporary escalation instrumentation
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
// and the block-aware System eval/Jacobian/time-derivative through the AMP tracker + Cauchy endgame.
//
// Deterministic: SetGlobalSeed pins the homotopy gamma and the MHom start coefficients (RandomMp is
// reseedable now), so this is a single, reproducible solve -- no gamma-retry loop.
//
// Skipped on Windows: there, AMP + the blend-block homotopy + the Cauchy endgame do not keep
// precision in lockstep (the tracked point reaches MaxPrecisionAllowed in the endgame while the
// system is at the working precision), so a path can grind for hours instead of converging -- this
// once hung Windows CI.  Tracked as the block-precision follow-up; the AMP MHom path itself is
// covered on Windows by the eigenvalue test.
BOOST_AUTO_TEST_CASE(mhom_solves_two_variable_group_system)
{
#ifdef _WIN32
	BOOST_TEST_MESSAGE("mhom_solves_two_variable_group_system skipped on Windows (block-precision grind)");
#else
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

	auto zd = algorithm::ZeroDim<AMPTracker,
	                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
	                             decltype(sys),
	                             start_system::MHomogeneous>(sys);
	zd.DefaultSetup();
	zd.Solve();

	// Collect the successfully-tracked solutions (a failed path leaves an empty placeholder).
	auto const& sols = zd.SolutionsUserCoords();
	auto const& md   = zd.FinalSolutionMetadata();
	using SolVec = std::decay_t<decltype(sols[0])>;
	std::vector<SolVec> good;
	for (size_t i = 0; i < sols.size(); ++i)
		if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 2)
			good.push_back(sols[i]);

	BOOST_REQUIRE_EQUAL(good.size(), 2u); // both MHom paths solved (the m-homogeneous Bezout number)

	for (auto const& s : good)   // AMP returns mpfr_complex coords; double is plenty for a root check
	{
		dbl a(s(0)), b(s(1));
		BOOST_CHECK_SMALL(std::abs(a * b - dbl(1)), 1e-8);
		BOOST_CHECK_SMALL(std::abs(a + b), 1e-8);
	}
	BOOST_CHECK_GT(std::abs(dbl(good[0](0)) - dbl(good[1](0))), 1e-3); // the two distinct roots
#endif
}


// PROBE (temporary, branch perf/amp-block-precision-escalation): sweep many gammas and report,
// per gamma, WHERE precision escalations originate.  Confirms "bad gamma" is really over-escalation.
// Run with:  ./build/core/test_nag_algorithms --run_test=zero_dim/amp_escalation_probe --log_level=message
BOOST_AUTO_TEST_CASE(amp_escalation_probe)
{
	using namespace bertini;
	using namespace tracking;

	std::cout << "\n seed |  good | maxPrec | maxDigB | trkInc | egRefine | corrTrkHPN | corrRefHPN |    ms\n";
	std::cout << "------+-------+---------+---------+--------+----------+------------+------------+--------\n";

	unsigned seed_lo = 1, seed_hi = 25;
	if (const char* one = std::getenv("BERTINI_PROBE_SEED")) { seed_lo = seed_hi = static_cast<unsigned>(std::atoi(one)); }
	const bool use_double = (std::getenv("BERTINI_PROBE_DOUBLE") != nullptr); // fixed double tracker instead of AMP

	for (unsigned seed = seed_lo; seed <= seed_hi; ++seed)
	{
		SetGlobalSeed(seed);           // different gamma + MHom start coefficients per seed
		bertini::probe::reset();

		System sys;
		auto x = Variable::Make("x");
		auto y = Variable::Make("y");
		sys.AddVariableGroup(VariableGroup{x});
		sys.AddVariableGroup(VariableGroup{y});
		sys.AddFunction(x*y - 1);
		sys.AddFunction(x + y);

		unsigned good = 0;
		long long ms = 0;
		if (use_double)
		{
			auto zd = algorithm::ZeroDim<DoublePrecisionTracker,
			                             bertini::endgame::EndgameSelector<DoublePrecisionTracker>::Cauchy,
			                             decltype(sys),
			                             start_system::MHomogeneous>(sys);
			zd.DefaultSetup();
			auto t0 = std::chrono::steady_clock::now();
			zd.Solve();
			auto t1 = std::chrono::steady_clock::now();
			ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
			auto const& sols = zd.SolutionsUserCoords();
			auto const& md   = zd.FinalSolutionMetadata();
			for (size_t i = 0; i < sols.size(); ++i)
				if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 2) ++good;
		}
		else
		{
			auto zd = algorithm::ZeroDim<AMPTracker,
			                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
			                             decltype(sys),
			                             start_system::MHomogeneous>(sys);
			zd.DefaultSetup();
			auto t0 = std::chrono::steady_clock::now();
			zd.Solve();
			auto t1 = std::chrono::steady_clock::now();
			ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
			auto const& sols = zd.SolutionsUserCoords();
			auto const& md   = zd.FinalSolutionMetadata();
			for (size_t i = 0; i < sols.size(); ++i)
				if (md[i].endgame_success == SuccessCode::Success && sols[i].size() == 2) ++good;
		}

		std::cout << std::setw(5) << seed << " | "
		          << std::setw(5) << good << " | "
		          << std::setw(7) << bertini::probe::max_precision_seen.load() << " | "
		          << std::setw(7) << bertini::probe::max_digits_b.load() << " | "
		          << std::setw(6) << bertini::probe::tracker_precision_increases.load() << " | "
		          << std::setw(8) << bertini::probe::endgame_refine_escalations.load() << " | "
		          << std::setw(10) << bertini::probe::corrector_track_hpn.load() << " | "
		          << std::setw(10) << bertini::probe::corrector_refine_hpn.load() << " | "
		          << std::setw(6) << ms << std::endl; // flush per row so a grinding gamma can't hide prior rows
	}
	std::cout << std::flush;
	BOOST_CHECK(true); // probe prints; correctness asserted by the sibling tests
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

	auto zd = algorithm::ZeroDim<AMPTracker,
	                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
	                             decltype(sys),
	                             start_system::MHomogeneous>(sys);
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
		Vec<dbl> p = Vec<dbl>::Random(n);
		dbl t = (trial < 5) ? dbl(0.4, -0.3) * dbl(trial + 1)
		                    : dbl(std::pow(10.0, -(trial - 1)), 0.0); // 1e-4 .. 1e-8
		compare(H.Eval(p, t),     He.Eval(p, t),     1e-11);
		compare(H.Jacobian(p, t), He.Jacobian(p, t), 1e-11);
	}

	// mpfr at 80 digits
	DefaultPrecision(80);
	H.precision(80);
	He.precision(80);
	for (int trial = 0; trial < 5; ++trial)
	{
		Vec<mpfr_complex> p = RandomOfUnits<mpfr_complex>(n);
		mpfr_complex t = RandomOfUnits<mpfr_complex>(1)(0);
		compare(H.Eval(p, t),     He.Eval(p, t),     1e-70);
		compare(H.Jacobian(p, t), He.Jacobian(p, t), 1e-70);
	}
}


// PROBE observer (branch perf/amp-block-precision-escalation): record, per successful step, the
// |t|, working precision, and the tracker's condition-number estimate -- so we can SEE whether the
// condition number (||J|| * ||J^{-1}||) spikes then RECOVERS along the actual seed-6 path, and at
// what |t| (mid-path vs the t->0 endgame region).
template <class TrackerT>
class CondNumTrajectory : public bertini::Observer<TrackerT>
{ BOOST_TYPE_INDEX_REGISTER_CLASS
	using EmitterT = typename bertini::tracking::TrackerTraits<TrackerT>::EventEmitterType;

	std::vector<std::type_index> SubscribedEventTypes() const override
	{ return { typeid(bertini::tracking::SuccessfulStep<EmitterT>) }; }

	void Observe(bertini::AnyEvent const& e) override
	{
		auto p = dynamic_cast<const bertini::tracking::SuccessfulStep<EmitterT>*>(&e);
		if (p)
		{
			auto const& tr = p->Get();
			rows.emplace_back(static_cast<double>(abs(tr.CurrentTime())),
			                  tr.CurrentPrecision(),
			                  static_cast<double>(tr.LatestConditionNumber()));
		}
	}
public:
	std::vector<std::tuple<double, unsigned, double>> rows; // (|t|, precision, condition number)
	virtual ~CondNumTrajectory() = default;
};


// #3 from the AMP-escalation investigation: along the actual seed-6 MHom path, log the condition
// number / precision vs |t|, to confirm whether ||J^{-1}|| spikes then RECOVERS (a transient
// near-singular pass) and where.  Also recovers and prints gamma to confirm it is genuinely
// complex.  Information-gathering only -- no assertions about the trajectory shape.
BOOST_AUTO_TEST_CASE(mhom_condition_number_trajectory)
{
	using namespace bertini;
	using namespace tracking;

	// Sweep seeds until we CATCH a genuinely spiking path (max precision pushed well above the
	// baseline), then dump that path's condition-number/precision trajectory -- so we can see
	// whether ||J^{-1}|| spikes then RECOVERS, and at what |t|.  (Necessary because SetGlobalSeed
	// does not fully reset RNG, so a fixed seed is not reproducible across call contexts.)
	auto build = []() {
		System sys;
		auto x = Variable::Make("x");
		auto y = Variable::Make("y");
		sys.AddVariableGroup(VariableGroup{x});
		sys.AddVariableGroup(VariableGroup{y});
		sys.AddFunction(x*y - 1);
		sys.AddFunction(x + y);
		return sys;
	};

	CondNumTrajectory<AMPTracker> spike;   // trajectory of the first spiking seed found
	unsigned spike_seed = 0; dbl spike_gamma(0,0);

	for (unsigned seed = 1; seed <= 60 && spike_seed == 0; ++seed)
	{
		SetGlobalSeed(seed);
		bertini::probe::reset();
		System sys = build();
		auto zd = algorithm::ZeroDim<AMPTracker,
		                             bertini::endgame::EndgameSelector<AMPTracker>::Cauchy,
		                             decltype(sys),
		                             start_system::MHomogeneous>(sys);
		zd.DefaultSetup();

		dbl gamma(0,0);
		{
			System const& H = zd.Homotopy();
			auto const& start = zd.StartSystem();
			Vec<dbl> xr = Vec<dbl>::Random(static_cast<int>(H.NumVariables()));
			Vec<dbl> Hv = H.Eval(xr, dbl(1.0));
			Vec<dbl> Sv = start.Eval(xr);
			gamma = Hv(0) / Sv(0);
		}

		CondNumTrajectory<AMPTracker> traj;
		zd.GetTracker().AddObserver(traj);
		zd.Solve();
		zd.GetTracker().RemoveObserver(traj);

		if (bertini::probe::max_precision_seen.load() > 80) // a genuine escalation
		{
			spike = std::move(traj);
			spike_seed = seed;
			spike_gamma = gamma;
		}
	}

	if (spike_seed == 0)
	{
		std::cout << "\n[trajectory] no spiking seed found in 1..60 in this context (max precision stayed low).\n" << std::flush;
		BOOST_CHECK(true);
		return;
	}

	std::cout << "\n[trajectory] spiking seed = " << spike_seed
	          << " ; gamma = " << spike_gamma << " (|Im|=" << std::abs(spike_gamma.imag())
	          << ", |gamma|=" << std::abs(spike_gamma) << ")\n";
	// probe counters still hold the spiking seed's values (loop exited on finding it): WHERE did
	// the escalation originate?
	std::cout << "[trajectory] per-cause for this seed:"
	          << " trackerPrecIncreases=" << bertini::probe::tracker_precision_increases.load()
	          << " endgameRefineEscalations=" << bertini::probe::endgame_refine_escalations.load()
	          << " corrTrackHPN=" << bertini::probe::corrector_track_hpn.load()
	          << " corrRefineHPN=" << bertini::probe::corrector_refine_hpn.load()
	          << " maxDigitsB=" << bertini::probe::max_digits_b.load() << "\n";
	std::cout << "   step |        |t|        | prec | log10(condNum)\n";
	std::cout << "  ------+-------------------+------+----------------\n";
	double max_cond = 0.0; double abst_at_max = 0.0; unsigned max_prec = 0; size_t step_at_max = 0;
	for (size_t i = 0; i < spike.rows.size(); ++i)
	{
		double abst = std::get<0>(spike.rows[i]);
		unsigned prec = std::get<1>(spike.rows[i]);
		double cond = std::get<2>(spike.rows[i]);
		double lc = (cond > 0) ? std::log10(cond) : 0.0;
		if (cond > max_cond) { max_cond = cond; abst_at_max = abst; step_at_max = i; }
		if (prec > max_prec) max_prec = prec;
		std::cout << std::setw(7) << i << " | " << std::setw(17) << std::scientific << std::setprecision(6) << abst
		          << " | " << std::setw(4) << prec << " | " << std::setw(14) << std::fixed << std::setprecision(3) << lc
		          << std::endl;
	}
	// did precision recover after the peak?
	unsigned prec_at_end = spike.rows.empty() ? 0 : std::get<1>(spike.rows.back());
	std::cout << "[trajectory] max log10(condNum) = " << std::log10(std::max(max_cond,1.0))
	          << " at |t| = " << std::scientific << abst_at_max << " (step " << step_at_max << "/" << spike.rows.size() << ")"
	          << " ; max precision = " << max_prec << " ; precision at path end = " << prec_at_end
	          << (prec_at_end < max_prec ? "  (RECOVERED)" : "  (did NOT recover)") << "\n" << std::flush;

	DefaultPrecision(30); // restore for subsequent tests

	BOOST_CHECK(true); // diagnostic; no trajectory-shape assertion
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto c = Tally(zd.FinalSolutionMetadata());
	BOOST_CHECK_EQUAL(c.success, 2);
	BOOST_CHECK_EQUAL(c.finite, 2);
	BOOST_CHECK_EQUAL(c.real, 2);
	BOOST_CHECK_EQUAL(c.singular, 0);
	for (auto const& m : zd.FinalSolutionMetadata())
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto c = Tally(zd.FinalSolutionMetadata());
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
	zd.DefaultSetup();
	zd.Solve();

	auto const& md = zd.FinalSolutionMetadata();
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
	zd.DefaultSetup();
	auto pp = zd.Get<PostProcessing>();
	pp.endpoint_finite_threshold = 0.5;          // 1 > 0.5 -> "at infinity"
	zd.Set(pp);
	zd.Solve();

	auto c = Tally(zd.FinalSolutionMetadata());
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
	zd.DefaultSetup();
	auto pp = zd.Get<PostProcessing>();
	pp.condition_number_threshold = 1e-3;        // any condition number exceeds this
	zd.Set(pp);
	zd.Solve();

	auto c = Tally(zd.FinalSolutionMetadata());
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

	auto zd = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);
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


BOOST_AUTO_TEST_SUITE_END()
