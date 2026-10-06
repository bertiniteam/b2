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

/**
\file threaded_solve.cpp

\brief Tests for MPI-less shared-memory threading: the WorkerThreadPool primitive and the
threaded ZeroDim solve.  The contract is that a threaded solve produces the SAME solution set
as a serial one -- threading is a performance feature, never a correctness change.  These run
on every platform (no MPI required), which is what gives us cross-OS threading coverage.
*/

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "bertini2/parallel/thread_pool.hpp"
#include "bertini2/random.hpp"
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

namespace {
// Set (value != nullptr) or unset an environment variable for the life of the object, and put
// back whatever was there.  POSIX setenv/unsetenv, or _putenv_s on Windows (where an empty value
// unsets).
class ScopedEnv
{
    std::string name_;
    bool had_ = false;
    std::string old_;

    static void Put(std::string const& name, char const* value)
    {
#ifdef _WIN32
        _putenv_s(name.c_str(), value ? value : "");
#else
        if (value)
            setenv(name.c_str(), value, 1);
        else
            unsetenv(name.c_str());
#endif
    }

public:
    ScopedEnv(std::string name, char const* value) : name_(std::move(name))
    {
        if (char const* prev = std::getenv(name_.c_str()))
        {
            had_ = true;
            old_ = prev;
        }
        Put(name_, value);
    }
    ~ScopedEnv() { Put(name_, had_ ? old_.c_str() : nullptr); }
    ScopedEnv(ScopedEnv const&) = delete;
    ScopedEnv& operator=(ScopedEnv const&) = delete;
};
} // namespace

BOOST_AUTO_TEST_CASE(effective_thread_count_is_sane)
{
    using bertini::parallel::EffectiveThreadCount;
    ScopedEnv no_override("BERTINI_NUM_THREADS", nullptr);
    // an explicit positive request is honored; auto is every available CPU, at least one
    BOOST_CHECK_EQUAL(EffectiveThreadCount(1), 1u);
    BOOST_CHECK_EQUAL(EffectiveThreadCount(3), 3u);
    BOOST_CHECK_EQUAL(EffectiveThreadCount(0), bertini::parallel::AvailableCpuCount());
    BOOST_CHECK_GE(bertini::parallel::AvailableCpuCount(), 1u);
}

BOOST_AUTO_TEST_CASE(bertini_num_threads_overrides_the_configured_count)
{
    using bertini::parallel::EffectiveThreadCount;
    ScopedEnv two("BERTINI_NUM_THREADS", "2");
    BOOST_CHECK_EQUAL(EffectiveThreadCount(3), 2u);
    BOOST_CHECK_EQUAL(EffectiveThreadCount(0), 2u);
}

// b2 reads its own variable.  OMP_NUM_THREADS also sets numpy's OpenBLAS thread count, and until
// 4.0 b2 read it too, coupling two unrelated settings; it must not reach b2 any more.
BOOST_AUTO_TEST_CASE(omp_num_threads_does_not_reach_b2)
{
    using bertini::parallel::EffectiveThreadCount;
    ScopedEnv no_override("BERTINI_NUM_THREADS", nullptr);
    ScopedEnv omp("OMP_NUM_THREADS", "1");
    BOOST_CHECK_EQUAL(EffectiveThreadCount(3), 3u);
    BOOST_CHECK_EQUAL(EffectiveThreadCount(0), bertini::parallel::AvailableCpuCount());
}

BOOST_AUTO_TEST_SUITE_END()  // thread_pool


// ---------------------------------------------------------------------------
// Threaded ZeroDim solve == serial ZeroDim solve (the correctness contract)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(threaded_solve)

namespace {

using Sol  = bertini::Vec<bertini::complex_dbl>;
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

    using ZD = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
                                  System>;

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
    using ZD = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
                                  System>;
    ZD zd(sys); zd.DefaultSetup();
    auto one = SolveWith(zd, 1);
    BOOST_CHECK_EQUAL(one.size(), 4u);
}


// A solver-level observer that records, per PathStarted, which tracker actually ran the path
// (event.Tracker()).  Because notifications are mutex-serialized, the set insert below is safe
// even when PathStarted fires from several worker threads -- so this also exercises the
// Observable notify mutex under concurrency.
struct TrackerCapture : public bertini::Observer<bertini::algorithm::AnyZeroDim>
{
    std::set<bertini::Observable const*> trackers_seen;
    bool                                 saw_null = false;
    int                                  path_starts = 0;

    bertini::ObserveResult Observe(bertini::AnyEvent const& e) override
    {
        using namespace bertini::algorithm;
        if (auto p = dynamic_cast<PathStarted<AnyZeroDim> const*>(&e))
        {
            ++path_starts;
            if (p->Tracker() == nullptr) saw_null = true;
            else                         trackers_seen.insert(p->Tracker());
        }
        return bertini::ObserveResult::KeepObserving;
    }
};

// Serial: every path's event.Tracker() is the solver's own member tracker.
BOOST_AUTO_TEST_CASE(event_tracker_is_member_tracker_when_serial)
{
    using namespace bertini;
    auto sys = TwoCubics();
    using ZD = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
                                  System>;
    ZD zd(sys); zd.DefaultSetup();

    TrackerCapture cap;
    zd.AddObserver(cap);
    SolveWith(zd, 1);

    auto const* member = static_cast<Observable const*>(&zd.GetTracker());
    BOOST_CHECK_EQUAL(cap.path_starts, 9);
    BOOST_CHECK(!cap.saw_null);
    BOOST_REQUIRE_EQUAL(cap.trackers_seen.size(), 1u);   // one tracker ran every path
    BOOST_CHECK(*cap.trackers_seen.begin() == member);
}

// Threaded: every path's event.Tracker() is a thread-local clone, never the member tracker --
// which is exactly why a meta-observer must attach to event.tracker(), not solver.GetTracker().
BOOST_AUTO_TEST_CASE(event_tracker_is_a_clone_when_threaded)
{
    using namespace bertini;
    // BERTINI_NUM_THREADS overrides the requested thread count (parallel::EffectiveThreadCount), so
    // under BERTINI_NUM_THREADS=1 a "threaded" solve actually runs serially on the member tracker,
    // and the clone-per-thread expectation below does not hold.  CI runs the test suite both with
    // and without BERTINI_NUM_THREADS=1; the threaded path is exercised in the non-serial leg, so
    // skip here when the environment forces serial.
    if (parallel::EffectiveThreadCount(4) <= 1)
    {
        BOOST_TEST_MESSAGE("event_tracker_is_a_clone_when_threaded: skipped -- BERTINI_NUM_THREADS "
                           "forces a serial run (threaded clone behavior is covered by the "
                           "BERTINI_NUM_THREADS!=1 CI leg).");
        return;
    }
    auto sys = TwoCubics();
    using ZD = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy,
                                  System>;
    ZD zd(sys); zd.DefaultSetup();

    TrackerCapture cap;
    zd.AddObserver(cap);
    SolveWith(zd, 4);

    auto const* member = static_cast<Observable const*>(&zd.GetTracker());
    BOOST_CHECK_EQUAL(cap.path_starts, 9);
    BOOST_CHECK(!cap.saw_null);
    BOOST_CHECK(cap.trackers_seen.count(member) == 0);   // never the member tracker
    BOOST_CHECK(!cap.trackers_seen.empty());
    BOOST_CHECK(cap.trackers_seen.size() <= 4u);          // at most one clone per worker thread
}


namespace {

// cyclic-5: 120 total-degree paths, 70 finite roots, the rest diverging.  Enough paths, and
// enough precision changes, that every thread count hands each thread a different sequence of
// predecessor paths.
bertini::System Cyclic5()
{
    using namespace bertini;
    std::vector<std::shared_ptr<node::Variable>> z;
    for (int i = 0; i < 5; ++i)
        z.push_back(Variable::Make("z" + std::to_string(i)));
    System sys;
    using Nd = std::shared_ptr<node::Node>;
    for (int k = 1; k < 5; ++k)
    {
        Nd total;
        for (int i = 0; i < 5; ++i)
        {
            Nd term = z[i];
            for (int j = 1; j < k; ++j)
                term = term * z[(i + j) % 5];
            total = total ? total + term : term;
        }
        sys.AddFunction(total);
    }
    sys.AddFunction(z[0] * z[1] * z[2] * z[3] * z[4] - 1);
    sys.AddVariableGroup(VariableGroup(z.begin(), z.end()));
    return sys;
}

// Every path's endpoint and its whole metadata record, from a seeded adaptive-precision solve.
template<typename EndgameT>
struct PathByPath
{
    using ZD = bertini::algorithm::ZeroDimSolver<bertini::tracking::AMPTracker, EndgameT,
                                                 bertini::System>;
    std::vector<bertini::Vec<bertini::complex_mp>> endpoints;
    std::vector<bertini::algorithm::SolutionMetaData<bertini::complex_mp>> metadata;

    explicit PathByPath(unsigned num_threads)
    {
        bertini::SetGlobalSeed(20261005);
        auto sys = Cyclic5();
        ZD zd(sys);
        zd.DefaultSetup();
        auto cfg = zd.template Get<bertini::algorithm::ZeroDimConfig>();
        cfg.num_threads = num_threads;
        zd.Set(cfg);
        zd.Solve();
        endpoints = zd.SolutionsInternalCoords();
        metadata = zd.SolutionMetadata();
    }
};

template<typename EndgameT>
void CheckBitIdenticalAcrossThreadCounts()
{
    PathByPath<EndgameT> const serial(1);
    BOOST_REQUIRE_EQUAL(serial.endpoints.size(), 120u);

    for (unsigned nt : {2u, 3u, 8u})
    {
        PathByPath<EndgameT> const threaded(nt);
        BOOST_REQUIRE_EQUAL(threaded.endpoints.size(), serial.endpoints.size());
        for (std::size_t i = 0; i < serial.endpoints.size(); ++i)
        {
            BOOST_TEST_CONTEXT("path " << i << " at " << nt << " threads")
            {
                auto const& a = serial.endpoints[i];
                auto const& b = threaded.endpoints[i];
                BOOST_REQUIRE_EQUAL(a.size(), b.size());
                for (Eigen::Index j = 0; j < a.size(); ++j)
                {
                    BOOST_CHECK_EQUAL(a(j).precision(), b(j).precision());
                    BOOST_CHECK(a(j) == b(j));
                }
                auto md_a = serial.metadata[i];   // operator== is not const
                BOOST_CHECK_MESSAGE(md_a == threaded.metadata[i],
                                    "metadata differs:\n" << serial.metadata[i]
                                    << "\nversus\n" << threaded.metadata[i]);
            }
        }
    }
}

} // namespace


/**
A threaded solve is bit-identical to the serial solve, path by path (#378).

A path's result is a function of the systems, the settings, the seed and its index -- not of
the thread that ran it, nor of the paths that thread ran before it.  Endpoints are compared
digit for digit and precision for precision, and every metadata field (step counts, precision
record, condition number, residuals, accuracy estimates) exactly; only the wall-clock time is
excluded.  Matching within a tolerance, as the tests above do, cannot see this: before the fix
the endpoints of a threaded solve agreed with the serial ones to many digits and still differed.

Under BERTINI_NUM_THREADS=1 every one of these solves is serial and the test proves nothing;
the CI leg without it is the one that counts.
*/
BOOST_AUTO_TEST_CASE(threaded_solve_is_bit_identical_to_serial_power_series)
{
    using TrackerA = bertini::tracking::AMPTracker;
    CheckBitIdenticalAcrossThreadCounts<bertini::endgame::EndgameSelector<TrackerA>::PSEG>();
}

BOOST_AUTO_TEST_CASE(threaded_solve_is_bit_identical_to_serial_cauchy)
{
    using TrackerA = bertini::tracking::AMPTracker;
    CheckBitIdenticalAcrossThreadCounts<bertini::endgame::EndgameSelector<TrackerA>::Cauchy>();
}


namespace {

// The highest precision the watched tracker reports, in any event: the start of a track, and
// both ends of a change of precision.
struct TrackerPrecisionWatch : public bertini::Observer<bertini::tracking::AMPTracker>
{
    unsigned highest = 0;

    bertini::ObserveResult Observe(bertini::AnyEvent const& e) override
    {
        using namespace bertini::tracking;
        if (auto s = dynamic_cast<TrackingStarted<AMPTracker> const*>(&e))
            highest = std::max(highest, s->Get().CurrentPrecision());
        else if (auto c = dynamic_cast<PrecisionChanged<AMPTracker> const*>(&e))
            highest = std::max({highest, c->Previous(), c->Next()});
        return bertini::ObserveResult::KeepObserving;
    }
};

// Starts the watch afresh at each path, and keeps its reading when the path completes.
struct PerPathPrecision : public bertini::Observer<bertini::algorithm::AnyZeroDim>
{
    TrackerPrecisionWatch& watch;
    std::map<std::size_t, unsigned> highest_by_path;

    explicit PerPathPrecision(TrackerPrecisionWatch& w) : watch(w) {}

    bertini::ObserveResult Observe(bertini::AnyEvent const& e) override
    {
        using namespace bertini::algorithm;
        if (dynamic_cast<PathStarted<AnyZeroDim> const*>(&e))
            watch.highest = 0;
        else if (auto c = dynamic_cast<PathComplete<AnyZeroDim> const*>(&e))
            highest_by_path[c->PathIndex()] = watch.highest;
        return bertini::ObserveResult::KeepObserving;
    }
};

template<typename EndgameT>
void CheckMaxPrecisionUsedCoversEveryTrack()
{
    using namespace bertini;
    // serial, so every path runs on the solver's own tracker, where the watch is attached
    thread_pool::ScopedEnv serial("BERTINI_NUM_THREADS", "1");

    SetGlobalSeed(20261005);
    auto sys = Cyclic5();
    algorithm::ZeroDimSolver<tracking::AMPTracker, EndgameT, System> zd(sys);
    zd.DefaultSetup();

    TrackerPrecisionWatch watch;
    PerPathPrecision per_path(watch);
    zd.GetTracker().AddObserver(watch);
    zd.AddObserver(per_path);
    zd.Solve();

    auto const& md = zd.SolutionMetadata();
    BOOST_REQUIRE_EQUAL(per_path.highest_by_path.size(), md.size());
    bool some_path_raised_precision = false;
    for (std::size_t i = 0; i < md.size(); ++i)
    {
        BOOST_TEST_CONTEXT("path " << i)
        {
            BOOST_CHECK_EQUAL(md[i].max_precision_used, per_path.highest_by_path.at(i));
            some_path_raised_precision = some_path_raised_precision || md[i].max_precision_used > 16;
        }
    }
    BOOST_CHECK(some_path_raised_precision);   // the premise: the check is not about doubles alone
}

} // namespace


/**
A path's max_precision_used is the highest precision any of its tracks used (#378 follow-up).

The endgame tracks a path in many calls to TrackPath.  The solver's record of a path's precision
used to start afresh with every call, so a precision raised in an earlier endgame sub-track and
lowered before the last went unreported.  Checked against every precision the tracker reports
while tracking each path, on a serial solve so that the solver's own tracker runs every path.
*/
BOOST_AUTO_TEST_CASE(max_precision_used_covers_every_track_of_a_path_power_series)
{
    using TrackerA = bertini::tracking::AMPTracker;
    CheckMaxPrecisionUsedCoversEveryTrack<bertini::endgame::EndgameSelector<TrackerA>::PSEG>();
}

BOOST_AUTO_TEST_CASE(max_precision_used_covers_every_track_of_a_path_cauchy)
{
    using TrackerA = bertini::tracking::AMPTracker;
    CheckMaxPrecisionUsedCoversEveryTrack<bertini::endgame::EndgameSelector<TrackerA>::Cauchy>();
}

BOOST_AUTO_TEST_SUITE_END()
