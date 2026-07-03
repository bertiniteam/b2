//This file is part of Bertini 2.
//
//zero_dim_records.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//zero_dim_records.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with zero_dim_records.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file The records seam (ADR-0046): solve() is ensure-answered.  A recording solve emits
a run header + one track record per path; an identical ask against the same directory
hydrates instead of computing; a PARTIAL directory (the kill-and-rerun case) hydrates
what exists and computes only the rest; BERTINI_RECORDS_DIR attaches ambient records
with zero API calls.  Seed-rooted randomness (ADR-0044) is what makes the rebuilt
homotopy identical, so hydrated and computed results are directly comparable.
*/

#include <cstdlib>
#include <filesystem>
#include <fstream>

#include <boost/test/unit_test.hpp>

#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/records/output_directory.hpp"

using namespace bertini;
using Variable = node::Variable;
using TrackerT = tracking::DoublePrecisionTracker;
using ZD = algorithm::ZeroDimSolver<TrackerT, endgame::EndgameSelector<TrackerT>::Cauchy, System>;
namespace fs = std::filesystem;

namespace {

fs::path FreshDir(std::string const& name)
{
	auto path = fs::temp_directory_path() / ("b2_records_seam_" + name);
	fs::remove_all(path);
	return path;
}

// {x^2 - 1, y^2 - 1}: four well-separated roots; 4 total-degree paths
System TwoQuadrics()
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddFunction(pow(y, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});
	return sys;
}

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(zero_dim_records)

BOOST_AUTO_TEST_CASE(recording_solve_then_full_hydration)
{
	auto const dir = FreshDir("hydrate");

	// first solve: records everything
	SetGlobalSeed(42);
	auto sys_a = TwoQuadrics();
	ZD a(sys_a);
	a.DefaultSetup();
	a.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	a.Solve();
	BOOST_CHECK_EQUAL(a.NumPathsHydrated(), 0u);

	unsigned runs = 0, tracks = 0;
	for (auto const& rec : a.Records()->Scan())
	{
		auto const kind = std::string(rec.at("kind").as_string());
		if (kind == "run")
		{
			++runs;
			// the archive says which software wrote it (descriptive, not identity)
			auto const& producer = rec.at("producer").as_object();
			BOOST_CHECK_EQUAL(std::string(producer.at("name").as_string()), "bertini2");
			BOOST_CHECK(!producer.at("version").as_string().empty());
			BOOST_CHECK(!producer.at("commit").as_string().empty());
		}
		if (kind == "track")
		{
			++tracks;
			BOOST_CHECK_EQUAL(std::string(rec.at("status").as_string()), "success");
		}
	}
	BOOST_CHECK_EQUAL(runs, 1u);
	BOOST_CHECK_EQUAL(tracks, 4u);

	// identical ask (same seed => same homotopy, ADR-0044): hydrates, computes nothing
	SetGlobalSeed(42);
	auto sys_b = TwoQuadrics();
	ZD b(sys_b);
	b.DefaultSetup();
	b.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	b.Solve();
	BOOST_CHECK_EQUAL(b.NumPathsHydrated(), 4u);

	// hydrated state is identical to computed state (same installer, same records)
	auto const& first = a.SolutionsInternalCoords();
	auto const& again = b.SolutionsInternalCoords();
	BOOST_REQUIRE_EQUAL(first.size(), again.size());
	for (std::size_t ii = 0; ii < first.size(); ++ii)
		BOOST_CHECK_SMALL((first[ii] - again[ii]).norm(), 1e-14);
}

BOOST_AUTO_TEST_CASE(partial_directory_resumes_computing_only_the_missing)
{
	auto const full_dir = FreshDir("resume_full");
	auto const partial_dir = FreshDir("resume_partial");

	SetGlobalSeed(42);
	auto sys_full = TwoQuadrics();
	ZD full(sys_full);
	full.DefaultSetup();
	full.RecordTo(std::make_shared<records::OutputDirectory>(full_dir));
	full.Solve();

	// build the partial directory: the run header + only paths 0 and 2 (a crash after
	// two of four paths, in effect)
	{
		records::OutputDirectory partial(partial_dir);
		for (auto const& rec : full.Records()->Scan())
		{
			auto const kind = std::string(rec.at("kind").as_string());
			if (kind == "run")
				partial.Append(rec);
			else if (kind == "track")
			{
				auto const idx = rec.at("index").as_int64();
				if (idx == 0 || idx == 2)
					partial.Append(rec);
			}
		}
		// the target definition rides along too
		auto const header = [&full]{
			for (auto const& rec : full.Records()->Scan())
				if (std::string(rec.at("kind").as_string()) == "run")
					return rec;
			return boost::json::object{};
		}();
		auto const target_id = std::string(header.at("target_object").as_string());
		partial.PutDefinition(full.Records()->GetDefinition(target_id), target_id);
	}

	// the rerun: hydrates 2, computes 2, and the answers match the uninterrupted solve
	SetGlobalSeed(42);
	auto sys_resumed = TwoQuadrics();
	ZD resumed(sys_resumed);
	resumed.DefaultSetup();
	resumed.RecordTo(std::make_shared<records::OutputDirectory>(partial_dir));
	resumed.Solve();
	BOOST_CHECK_EQUAL(resumed.NumPathsHydrated(), 2u);

	auto const& expected = full.SolutionsInternalCoords();
	auto const& actual = resumed.SolutionsInternalCoords();
	BOOST_REQUIRE_EQUAL(expected.size(), actual.size());
	for (std::size_t ii = 0; ii < expected.size(); ++ii)
		BOOST_CHECK_SMALL((expected[ii] - actual[ii]).norm(), 1e-14);

	// and the partial directory is now complete: a further rerun hydrates everything
	SetGlobalSeed(42);
	auto sys_again = TwoQuadrics();
	ZD again(sys_again);
	again.DefaultSetup();
	again.RecordTo(std::make_shared<records::OutputDirectory>(partial_dir));
	again.Solve();
	BOOST_CHECK_EQUAL(again.NumPathsHydrated(), 4u);
}

BOOST_AUTO_TEST_CASE(different_seed_is_a_different_ask)
{
	auto const dir = FreshDir("seedsplit");

	SetGlobalSeed(42);
	auto sys_a = TwoQuadrics();
	ZD a(sys_a);
	a.DefaultSetup();
	a.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	a.Solve();

	SetGlobalSeed(43);   // different homotopy instance: nothing to hydrate
	auto sys_b = TwoQuadrics();
	ZD b(sys_b);
	b.DefaultSetup();
	b.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	b.Solve();
	BOOST_CHECK_EQUAL(b.NumPathsHydrated(), 0u);

	unsigned runs = 0;
	for (auto const& rec : b.Records()->Scan())
		if (std::string(rec.at("kind").as_string()) == "run")
			++runs;
	BOOST_CHECK_EQUAL(runs, 2u);   // two asks, two run headers, one directory
}

// Named regression (user, 2026-07-03): a path the endgame truncates near infinity is a
// VERDICT, not a failure -- the records must say "diverged", never "failed", or a
// cyclic5 audit shows 50 phantom failures.  And diverged paths hydrate like any other:
// every tracked path is in the store, so fails/divergences are auditable and memoized.
BOOST_AUTO_TEST_CASE(diverged_paths_are_recorded_as_diverged_not_failed)
{
	auto const dir = FreshDir("diverged");

	// {x*y - 1, x^2 - 1}: total degree 4 paths, exactly 2 finite solutions
	// ((1,1) and (-1,-1)) -- the other 2 paths head to infinity
	SetGlobalSeed(42);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System sys;
	sys.AddFunction(x * y - 1);
	sys.AddFunction(pow(x, 2) - 1);
	sys.AddVariableGroup(VariableGroup{x, y});

	ZD zd(sys);
	zd.DefaultSetup();
	zd.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	zd.Solve();

	unsigned successes = 0, diverged = 0, failed = 0;
	for (auto const& rec : zd.Records()->Scan())
	{
		if (std::string(rec.at("kind").as_string()) != "track")
			continue;
		auto const status = std::string(rec.at("status").as_string());
		if (status == "success") ++successes;
		if (status == "diverged") ++diverged;
		if (status == "failed") ++failed;
		// the durable rendering: code NAMES beside the integers
		BOOST_CHECK(rec.contains("endgame_success_code_name"));
		BOOST_CHECK(rec.contains("pre_endgame_success_code_name"));
	}
	BOOST_CHECK_EQUAL(successes, 2u);
	BOOST_CHECK_EQUAL(diverged, 2u);
	BOOST_CHECK_EQUAL(failed, 0u);

	// diverged paths are answers: a rerun hydrates ALL of them, recomputing none
	SetGlobalSeed(42);
	auto x2 = Variable::Make("x");
	auto y2 = Variable::Make("y");
	System sys_again;
	sys_again.AddFunction(x2 * y2 - 1);
	sys_again.AddFunction(pow(x2, 2) - 1);
	sys_again.AddVariableGroup(VariableGroup{x2, y2});
	ZD again(sys_again);
	again.DefaultSetup();
	again.RecordTo(std::make_shared<records::OutputDirectory>(dir));
	again.Solve();
	BOOST_CHECK_EQUAL(again.NumPathsHydrated(), 4u);
}

BOOST_AUTO_TEST_CASE(ambient_records_attach_from_the_environment)
{
	auto const dir = FreshDir("ambient");
#ifdef _WIN32
	_putenv_s("BERTINI_RECORDS_DIR", dir.string().c_str());
#else
	setenv("BERTINI_RECORDS_DIR", dir.string().c_str(), 1);
#endif

	SetGlobalSeed(7);
	auto sys = TwoQuadrics();
	ZD zd(sys);
	zd.DefaultSetup();
	zd.Solve();   // no RecordTo: the environment supplies the directory

#ifdef _WIN32
	_putenv_s("BERTINI_RECORDS_DIR", "");
#else
	unsetenv("BERTINI_RECORDS_DIR");
#endif

	BOOST_REQUIRE(zd.Records() != nullptr);
	BOOST_CHECK(fs::exists(dir / "README.txt"));
	BOOST_CHECK_EQUAL(zd.Records()->Scan().size(), 6u);   // 1 run + 4 tracks + 1 auto-declared result
}

BOOST_AUTO_TEST_SUITE_END()
