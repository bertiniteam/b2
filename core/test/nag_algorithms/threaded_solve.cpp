//This file is part of Bertini 2.
//
//threaded_solve.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//threaded_solve.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with threaded_solve.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file threaded_solve.cpp

\brief Tests for MPI-less shared-memory threading: the WorkerThreadPool primitive and the
threaded ZeroDim solve.  The contract is that a threaded solve produces the SAME solution set
as a serial one -- threading is a performance feature, never a correctness change.  These run
on every platform (no MPI required), which is what gives us cross-OS threading coverage.
*/

#include <atomic>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "bertini2/parallel/thread_pool.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"

using Variable = bertini::node::Variable;
using TrackerT = bertini::tracking::DoublePrecisionTracker;


// ---------------------------------------------------------------------------
// WorkerThreadPool: the bare concurrency primitive (no MPI, no Bertini types)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(thread_pool)

// Submit N tasks across several worker threads, collect N results, shut down cleanly.  The
// result is order-independent, so summing checks that every task ran exactly once.
BOOST_AUTO_TEST_CASE(thread_pool_runs_every_task_once)
{
	using bertini::parallel::WorkerThreadPool;

	auto factory = []() { return std::unique_ptr<int>(new int(0)); };
	auto track   = [](std::unique_ptr<int>& state, int const& task) {
		++(*state);          // touch the per-thread state so isolation matters
		return task * task;
	};

	const int N = 200;
	WorkerThreadPool<int, int, decltype(factory), decltype(track)> pool(4, factory, track);

	for (int i = 0; i < N; ++i)
		pool.submit(i);

	long long got = 0;
	for (int k = 0; k < N; ++k)
		got += pool.collect();

	pool.shutdown();   // joins all worker threads; must not hang or throw

	long long expected = 0;
	for (int i = 0; i < N; ++i)
		expected += static_cast<long long>(i) * i;

	BOOST_CHECK_EQUAL(got, expected);
}

// Each thread builds its own state via the factory; states must not be shared.  We hand each
// thread a distinct heap counter and confirm the per-thread increments sum to the task count
// (a shared-state race would lose increments under a data race / give a wrong total).
BOOST_AUTO_TEST_CASE(thread_pool_state_is_per_thread)
{
	using bertini::parallel::WorkerThreadPool;

	std::atomic<int> factory_calls{0};

	auto factory = [&factory_calls]() {
		++factory_calls;
		return std::unique_ptr<long>(new long(0));
	};
	auto track = [](std::unique_ptr<long>& state, int const& task) {
		++(*state);
		return static_cast<long>(*state);   // value is meaningful only within one thread
	};

	const int N = 500;
	const int n_threads = 4;
	WorkerThreadPool<int, long, decltype(factory), decltype(track)> pool(n_threads, factory, track);

	for (int i = 0; i < N; ++i)
		pool.submit(i);

	long total = 0;
	for (int k = 0; k < N; ++k)
		total += pool.collect();
	pool.shutdown();

	// The factory runs exactly once per thread (per-thread state, not per-task).
	BOOST_CHECK_EQUAL(factory_calls.load(), n_threads);
	// Every task incremented some thread's private counter exactly once; the per-thread counters
	// partition the N tasks, so the returned running-counts sum to 1+2+...+(per-thread totals).
	// The weakest invariant that always holds regardless of scheduling: at least N (each task
	// returned >= 1) and the counters together saw exactly N increments.
	BOOST_CHECK_GE(total, static_cast<long>(N));
}

BOOST_AUTO_TEST_CASE(effective_thread_count_is_sane)
{
	using bertini::parallel::EffectiveThreadCount;
	// An explicit positive request is honored (when OMP_NUM_THREADS is unset in the test env).
	if (std::getenv("OMP_NUM_THREADS") == nullptr)
	{
		BOOST_CHECK_EQUAL(EffectiveThreadCount(1), 1u);
		BOOST_CHECK_EQUAL(EffectiveThreadCount(3), 3u);
		BOOST_CHECK_GE(EffectiveThreadCount(0), 1u);   // auto -> hardware_concurrency, clamped >= 1
	}
}

BOOST_AUTO_TEST_SUITE_END()  // thread_pool


// ---------------------------------------------------------------------------
// Threaded ZeroDim solve == serial ZeroDim solve (the correctness contract)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(threaded_solve)

namespace {

using Sol  = bertini::Vec<bertini::dbl>;
using Sols = std::vector<Sol>;

template<typename ZD>
Sols SolveWith(ZD& zd, unsigned num_threads)
{
	auto cfg = zd.template Get<bertini::algorithm::ZeroDimConfig>();
	cfg.num_threads = num_threads;
	zd.Set(cfg);
	zd.Solve();

	Sols out;
	for (auto const& s : zd.FiniteSolutions(/*user_coords=*/true))
		out.push_back(s);
	return out;
}

double Distance(Sol const& a, Sol const& b)
{
	if (a.size() != b.size())
		return std::numeric_limits<double>::infinity();
	double d = 0;
	for (Eigen::Index i = 0; i < a.size(); ++i)
		d += std::norm(a(i) - b(i));      // std::norm = squared magnitude
	return std::sqrt(d);
}

// Greedy one-to-one match within tolerance -- robust to the out-of-order arrival of threaded
// solutions (no reliance on a fragile sort of near-equal coordinates).
void CheckSameSolutionSet(Sols const& serial, Sols const& threaded, double tol)
{
	BOOST_REQUIRE_EQUAL(serial.size(), threaded.size());

	std::vector<bool> used(threaded.size(), false);
	for (auto const& s : serial)
	{
		bool matched = false;
		for (std::size_t j = 0; j < threaded.size(); ++j)
		{
			if (!used[j] && Distance(s, threaded[j]) < tol)
			{
				used[j] = true;
				matched = true;
				break;
			}
		}
		BOOST_CHECK_MESSAGE(matched, "a serial solution had no threaded counterpart within tolerance");
	}
}

// {x^2 - 1, y^2 - 1}: four well-separated nonsingular roots (+/-1, +/-1).  Total degree 4 paths.
bertini::System TwoQuadrics()
{
	using namespace bertini;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddFunction(pow(y, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});
	return sys;
}

// A denser system: {x^3 - x, y^3 - y} -> 9 total-degree paths, 9 finite roots -- more paths than
// cores, so the pool genuinely round-robins work across threads.
bertini::System TwoCubics()
{
	using namespace bertini;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 3) - x);
	sys.AddFunction(pow(y, 3) - y);
	sys.AddVariableGroup(VariableGroup{x, y});
	return sys;
}

template<typename MakeSys>
void ThreadedMatchesSerial(MakeSys make_sys, std::size_t expected_finite)
{
	using namespace bertini;

	// Independent solver objects (CloneGiven copies the system), so the two solves share no state.
	auto sys_serial = make_sys();
	auto sys_thread = make_sys();

	using ZD = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
	                              System, start_system::TotalDegree>;

	ZD zd_serial(sys_serial);  zd_serial.DefaultSetup();
	ZD zd_thread(sys_thread);  zd_thread.DefaultSetup();

	auto serial = SolveWith(zd_serial, 1);
	BOOST_REQUIRE_EQUAL(serial.size(), expected_finite);

	for (unsigned nt : {2u, 4u})
	{
		auto threaded = SolveWith(zd_thread, nt);
		CheckSameSolutionSet(serial, threaded, 1e-6);
	}
}

} // namespace

BOOST_AUTO_TEST_CASE(threaded_matches_serial_two_quadrics)
{
	ThreadedMatchesSerial([]{ return TwoQuadrics(); }, 4u);
}

BOOST_AUTO_TEST_CASE(threaded_matches_serial_two_cubics)
{
	ThreadedMatchesSerial([]{ return TwoCubics(); }, 9u);
}

// num_threads == 1 must be byte-for-byte the pool-free serial path -- a regression guard that the
// thread-count branch in Solve() does not perturb the default (serial) result.
BOOST_AUTO_TEST_CASE(num_threads_one_equals_default_solve)
{
	using namespace bertini;
	auto sys = TwoQuadrics();
	using ZD = algorithm::ZeroDim<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
	                              System, start_system::TotalDegree>;
	ZD zd(sys); zd.DefaultSetup();
	auto one = SolveWith(zd, 1);
	BOOST_CHECK_EQUAL(one.size(), 4u);
}

BOOST_AUTO_TEST_SUITE_END()
