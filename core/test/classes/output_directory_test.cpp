//This file is part of Bertini 2.
//
//output_directory_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//output_directory_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with output_directory_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for the C++ structured output directory (b2rec/1, ADR-0045): ports of
the Python pilot's ledger tests, plus the cross-implementation bridge — this suite
WRITES an example directory (into the ctest working directory) that the Python
prototype's cross-impl test reads, and READS one the prototype writes when present.
*/

#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <thread>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/json.hpp>

#include "bertini2/records/output_directory.hpp"
#include "bertini2/detail/sha256.hpp"

using bertini::records::OutputDirectory;
namespace json = boost::json;
namespace fs = std::filesystem;

namespace {

fs::path FreshDir(std::string const& name)
{
	auto path = fs::temp_directory_path() / ("b2_outdir_test_" + name);
	fs::remove_all(path);
	return path;
}

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(output_directory)

BOOST_AUTO_TEST_CASE(construction_writes_the_self_documenting_readme)
{
	auto dir = FreshDir("readme");
	OutputDirectory out(dir);
	BOOST_CHECK(fs::exists(dir / "README.txt"));
	std::ifstream readme(dir / "README.txt");
	std::string text((std::istreambuf_iterator<char>(readme)), std::istreambuf_iterator<char>());
	BOOST_CHECK(text.find("b2rec/1") != std::string::npos);
	BOOST_CHECK(text.find("free to delete") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(definitions_are_idempotent_and_self_verifying)
{
	auto const dir = FreshDir("defs");
	OutputDirectory out(dir);
	auto const id1 = out.PutDefinition("hello records", "givens");
	auto const id2 = out.PutDefinition("hello records", "givens");
	BOOST_CHECK_EQUAL(id1, id2);
	BOOST_CHECK_EQUAL(out.GetDefinition(id1), "hello records");
	BOOST_CHECK_EQUAL(id1, bertini::detail::Sha256("hello records").Hex());

	auto const ext = out.PutDefinition("a rendering", "systems", std::string(64, 'a'));
	BOOST_CHECK_EQUAL(ext, std::string(64, 'a'));
	BOOST_CHECK(out.HasDefinition(ext));

	// layout: kind folder for the human, two-hex shard for scale, and the filename
	// carries kind + FULL digest + an honest extension (no concatenation to recover
	// the digest, and the file says what it is even away from its folder)
	BOOST_CHECK(fs::exists(dir / "definitions" / "givens" / id1.substr(0, 2)
	                       / ("given-" + id1 + ".txt")));
	BOOST_CHECK(fs::exists(dir / "definitions" / "systems" / "aa"
	                       / ("system-" + std::string(64, 'a') + ".txt")));
	auto const jid = out.PutDefinition("{\"kind\":\"probe\"}", "givens");
	BOOST_CHECK(fs::exists(dir / "definitions" / "givens" / jid.substr(0, 2)
	                       / ("given-" + jid + ".json")));
	// a role label weaves into the filename (presentation only; same resolution)
	auto const lid = out.PutDefinition("my input file", "givens", std::nullopt, "cli_input");
	BOOST_CHECK(fs::exists(dir / "definitions" / "givens" / lid.substr(0, 2)
	                       / ("given-cli_input-" + lid + ".txt")));
	BOOST_CHECK_EQUAL(out.GetDefinition(lid), "my input file");
	// ids resolve without the kind: readers follow bare ids from the records
	BOOST_CHECK_EQUAL(out.GetDefinition(ext), "a rendering");
}

BOOST_AUTO_TEST_CASE(append_scan_round_trips_in_one_date_named_session_file)
{
	auto dir = FreshDir("journal");
	OutputDirectory out(dir);
	out.Append({{"kind", "run"}, {"run", "runA"}, {"num_paths", 1}});
	out.Append({{"kind", "track"}, {"run", "runA"}, {"index", 0}, {"status", "success"}});

	auto const records = out.Scan();
	BOOST_REQUIRE_EQUAL(records.size(), 2u);
	BOOST_CHECK_EQUAL(std::string(records[0].at("kind").as_string()), "run");
	BOOST_CHECK_EQUAL(std::string(records[1].at("kind").as_string()), "track");

	unsigned journal_count = 0;
	std::string name;
	for (auto const& entry : fs::directory_iterator(dir / "history"))
		if (entry.path().extension() == ".jsonl")
		{
			++journal_count;
			name = entry.path().filename().string();
		}
	BOOST_CHECK_EQUAL(journal_count, 1u);
	BOOST_CHECK(name.size() > 8 && std::isdigit(name[0]));   // YYYYMMDD prefix
}

BOOST_AUTO_TEST_CASE(torn_final_line_is_tolerated_torn_interior_is_not)
{
	auto dir = FreshDir("torn");
	OutputDirectory out(dir);
	out.Append({{"kind", "run"}, {"run", "runB"}});
	// simulate a kill mid-append
	fs::path journal;
	for (auto const& entry : fs::directory_iterator(dir / "history"))
		if (entry.path().extension() == ".jsonl")
			journal = entry.path();
	{
		std::ofstream torn(journal, std::ios::app);
		torn << "{\"kind\": \"track\", \"ind";
	}
	BOOST_CHECK_EQUAL(out.Scan().size(), 1u);   // torn tail invisible to replay

	// a torn INTERIOR line is real corruption
	auto dir2 = FreshDir("torn_interior");
	OutputDirectory out2(dir2);
	{
		std::ofstream bad(dir2 / "history" / "20200101_000000-pid1.jsonl");
		bad << "{\"kind\": \"run\"\n{\"kind\": \"track\", \"index\": 0}\n";
	}
	BOOST_CHECK_THROW(out2.Scan(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(unknown_record_kinds_are_preserved)
{
	OutputDirectory out(FreshDir("openvocab"));
	out.Append({{"kind", "regen_start"}, {"level", 3}, {"whatever", "future"}});
	auto const records = out.Scan();
	BOOST_REQUIRE_EQUAL(records.size(), 1u);
	BOOST_CHECK_EQUAL(std::string(records[0].at("kind").as_string()), "regen_start");
	BOOST_CHECK_EQUAL(records[0].at("level").as_int64(), 3);
}

// The three truth stores separate concerns: history says what was asked, results/
// holds the computed paths (one file per run, header line first), definitions/ holds
// the inputs.  This exercises the results store's whole contract.
BOOST_AUTO_TEST_CASE(results_files_are_per_run_self_describing_and_torn_tolerant)
{
	auto const dir = FreshDir("results_store");
	OutputDirectory out(dir);

	out.EnsureResultsFile("abc123", {{"kind", "results_header"}, {"run", "abc123"}});
	out.EnsureResultsFile("abc123", {{"kind", "results_header"}, {"run", "DUPLICATE"}});  // idempotent
	out.AppendResult("abc123", {{"kind", "path"}, {"run", "abc123"}, {"index", 0},
	                            {"status", "success"}});
	out.AppendResult("abc123", {{"kind", "path"}, {"run", "abc123"}, {"index", 1},
	                            {"status", "diverged"}});

	// sharded like definitions/: results/<2 hex>/<run id>.jsonl
	auto const path = dir / "results" / "ab" / "abc123.jsonl";
	BOOST_REQUIRE(fs::exists(path));

	auto const records = out.ResultsOf("abc123");
	BOOST_REQUIRE_EQUAL(records.size(), 3u);
	BOOST_CHECK_EQUAL(std::string(records[0].at("kind").as_string()), "results_header");
	BOOST_CHECK_EQUAL(std::string(records[0].at("run").as_string()), "abc123");  // first header won
	BOOST_CHECK_EQUAL(records[1].at("index").as_int64(), 0);
	BOOST_CHECK_EQUAL(std::string(records[2].at("status").as_string()), "diverged");

	// a torn final line (kill mid-append) is invisible to replay
	{
		std::ofstream torn(path, std::ios::app);
		torn << "{\"kind\": \"path\", \"ind";
	}
	BOOST_CHECK_EQUAL(out.ResultsOf("abc123").size(), 3u);

	// nothing recorded for a run reads as empty, not an error
	BOOST_CHECK(out.ResultsOf("beef00").empty());

	// hostile run ids are refused outright: they become filesystem paths
	BOOST_CHECK_THROW(out.AppendResult("../../evil", {{"kind", "path"}}), std::invalid_argument);
	BOOST_CHECK_THROW(out.EnsureResultsFile("ABC123", {}), std::invalid_argument);  // not lowercase hex
	BOOST_CHECK(out.ResultsOf("../../../etc/passwd").empty());
}

BOOST_AUTO_TEST_CASE(index_renders_from_history_and_results)
{
	auto dir = FreshDir("views");
	OutputDirectory out(dir);
	auto const target_id = out.PutDefinition(
		"function f;\nvariable_group x, y;\nf = x^2+4*y^2-4;\n", "systems");

	out.EnsureResultsFile("abc123", {{"kind", "results_header"}, {"run", "abc123"}});
	out.AppendResult("abc123", {{"kind", "path"}, {"run", "abc123"}, {"index", 0},
	            {"status", "success"},
	            {"endpoint", json::array{json::array{"1.5", "0"}, json::array{"-0.66", "0"}}}});
	out.Append({{"kind", "run"}, {"schema", "b2rec/1"}, {"when", "2026-07-03 10:00"},
	            {"run", "abc123"}, {"op", "solve"}, {"target_digest", target_id},
	            {"target_rendering", "function f;\nvariable_group x, y;\nf = x^2+4*y^2-4;\n"},
	            {"num_paths", 1}});

	std::ifstream index(dir / "INDEX.txt");
	std::string index_text((std::istreambuf_iterator<char>(index)), std::istreambuf_iterator<char>());
	BOOST_CHECK(index_text.find("x^2+4*y^2-4") != std::string::npos);
	BOOST_CHECK(index_text.find("abc123") != std::string::npos);
	BOOST_CHECK(index_text.find("1/1 paths done") != std::string::npos);  // counted from results/

	// the retired monolith views never come back
	BOOST_CHECK(!fs::exists(dir / "results.json"));
	BOOST_CHECK(!fs::exists(dir / "RESULTS.txt"));
}

BOOST_AUTO_TEST_CASE(shared_is_one_instance_per_path_per_process)
{
	// a sweep attaching the ambient records to thousands of solvers must share ONE
	// writer (one session history file), not exhaust the per-second claim namespace
	auto dir = FreshDir("shared");
	auto a = OutputDirectory::Shared(dir);
	auto b = OutputDirectory::Shared(dir);
	BOOST_CHECK_EQUAL(a.get(), b.get());

	auto other = OutputDirectory::Shared(FreshDir("shared_other"));
	BOOST_CHECK(a.get() != other.get());

	// many attaches, one session file
	for (int i = 0; i < 100; ++i)
		OutputDirectory::Shared(dir)->Append({{"kind", "probe"}, {"i", i}});
	std::size_t session_files = 0;
	for (auto const& entry : fs::directory_iterator(dir / "history"))
		if (entry.path().extension() == ".jsonl")
			++session_files;
	BOOST_CHECK_EQUAL(session_files, 1u);

	// held weakly: once nobody references it, the instance is released
	std::weak_ptr<OutputDirectory> watch = a;
	a.reset(); b.reset();
	BOOST_CHECK(watch.expired());
}

BOOST_AUTO_TEST_CASE(hostile_ids_kinds_and_labels_are_refused)
{
	// these strings become filesystem paths: a hostile or buggy caller must not be
	// able to write outside definitions/ or smuggle separators into filenames
	auto const dir = FreshDir("hostile");
	OutputDirectory out(dir);
	BOOST_CHECK_THROW(out.PutDefinition("x", "systems", "../../evil"), std::invalid_argument);
	BOOST_CHECK_THROW(out.PutDefinition("x", "systems", std::string(64, 'Z')), std::invalid_argument);
	BOOST_CHECK_THROW(out.PutDefinition("x", "../escape"), std::invalid_argument);
	BOOST_CHECK_THROW(out.PutDefinition("x", "systems", std::nullopt, "a/b"), std::invalid_argument);
	BOOST_CHECK_THROW(out.PutDefinition("x", ""), std::invalid_argument);
	// and a hostile directory's RECORDS cannot point resolution outside definitions/
	BOOST_CHECK(!out.HasDefinition("../../../etc/passwd"));
	BOOST_CHECK(!out.HasDefinition(".."));
	BOOST_CHECK(!out.HasDefinition(""));
	BOOST_CHECK_THROW(out.GetDefinition("../../../etc/passwd"), std::runtime_error);
	// nothing escaped
	BOOST_CHECK(!fs::exists(dir.parent_path() / "evil"));
}

BOOST_AUTO_TEST_CASE(doubles_render_consistently_everywhere)
{
	// annotation values, results files, and history lines all use the same shortest
	// round-trip decimal form as the config views -- never uppercase-scientific
	auto dir = FreshDir("doubles");
	OutputDirectory out(dir);
	out.Append({{"kind", "run"}, {"run", "aa11"}, {"num_paths", 1}});
	out.AppendResult("aa11", {{"kind", "path"}, {"run", "aa11"}, {"index", 0},
	            {"status", "success"}, {"residual", 1e-05},
	            {"endpoint", json::array{json::array{"0.5", "0"}}}});
	out.Annotate("aa11", 0, "projection", json::value(0.5));
	out.Annotate("aa11", 0, "tolerance", json::value(1e-05));

	auto read = [](fs::path const& p) {
		std::ifstream in(p);
		return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	};
	auto const results_text = read(dir / "results" / "aa" / "aa11.jsonl");
	BOOST_CHECK(results_text.find("0.5") != std::string::npos);
	BOOST_CHECK(results_text.find("1e-05") != std::string::npos);
	BOOST_CHECK(results_text.find("5E-1") == std::string::npos);
	for (auto const& entry : fs::directory_iterator(dir / "history"))
		if (entry.path().extension() == ".jsonl")
		{
			auto const line_text = read(entry.path());
			BOOST_CHECK(line_text.find("5E-1") == std::string::npos);
			if (line_text.find("projection") != std::string::npos)
				BOOST_CHECK(line_text.find("0.5") != std::string::npos);
		}
	// round trip: the values come back as the same doubles
	for (auto const& rec : out.Scan())
		if (rec.if_contains("key") && rec.at("key").as_string() == "tolerance")
			BOOST_CHECK_EQUAL(rec.at("value").as_double(), 1e-05);
	for (auto const& rec : out.ResultsOf("aa11"))
		if (rec.if_contains("residual"))
			BOOST_CHECK_EQUAL(rec.at("residual").as_double(), 1e-05);
}

BOOST_AUTO_TEST_CASE(concurrent_appends_from_many_threads_all_land)
{
	// a Shared() instance may be written from several threads (ambient recording in
	// a threaded sweep): every record lands, none torn, one session file
	auto dir = FreshDir("threads");
	auto out = OutputDirectory::Shared(dir);
	constexpr int kThreads = 8, kEach = 50;
	std::vector<std::thread> workers;
	for (int t = 0; t < kThreads; ++t)
		workers.emplace_back([&out, t] {
			for (int i = 0; i < kEach; ++i)
				out->Append({{"kind", "probe"}, {"thread", t}, {"i", i}});
		});
	for (auto& w : workers)
		w.join();
	BOOST_CHECK_EQUAL(out->Scan().size(), static_cast<std::size_t>(kThreads * kEach));
	std::size_t session_files = 0;
	for (auto const& entry : fs::directory_iterator(dir / "history"))
		if (entry.path().extension() == ".jsonl")
			++session_files;
	BOOST_CHECK_EQUAL(session_files, 1u);
}

BOOST_AUTO_TEST_CASE(views_survive_dangling_and_garbage_references)
{
	// a run whose definition is MISSING, a result pointing at a nonexistent track,
	// and an annotation on a point that was never recorded: views must render
	// something, never throw, never crash
	auto dir = FreshDir("dangling");
	OutputDirectory out(dir);
	out.Append({{"kind", "run"}, {"run", "r1"}, {"num_paths", 1},
	            {"target_digest", std::string(64, '0')}});   // no such definition
	out.Append({{"kind", "result"}, {"name", "ghost"}, {"when", "now"},
	            {"points", json::array{json::object{{"run", "nope"}, {"index", 99}}}}});
	out.Annotate("never_ran", 7, "note", json::value("orphan"));
	BOOST_CHECK_NO_THROW(out.RefreshIndex());
	std::ifstream in(dir / "INDEX.txt");
	std::string index_text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	BOOST_CHECK(index_text.find("(unknown target)") != std::string::npos);
}

// The ONE place the BERTINI_RECORDS_DIR semantics live: unset = records on by
// default at ./bertini_output; a value = that directory; `none` (portable) or the
// empty string (POSIX-only: Windows deletes empty-valued variables) = off.
BOOST_AUTO_TEST_CASE(ambient_records_path_resolution)
{
	auto set_env = [](char const* value) {
#ifdef _WIN32
		_putenv_s("BERTINI_RECORDS_DIR", value);
#else
		setenv("BERTINI_RECORDS_DIR", value, 1);
#endif
	};

#ifdef _WIN32
	_putenv("BERTINI_RECORDS_DIR=");   // assigning empty DELETES the variable on Windows
#else
	unsetenv("BERTINI_RECORDS_DIR");
#endif
	auto const unset = bertini::records::AmbientRecordsPath();
	BOOST_REQUIRE(unset.has_value());
	BOOST_CHECK_EQUAL(*unset, "bertini_output");

	set_env("my_records");
	auto const chosen = bertini::records::AmbientRecordsPath();
	BOOST_REQUIRE(chosen.has_value());
	BOOST_CHECK_EQUAL(*chosen, "my_records");

	set_env("none");
	BOOST_CHECK(!bertini::records::AmbientRecordsPath().has_value());

#ifndef _WIN32
	setenv("BERTINI_RECORDS_DIR", "", 1);   // the POSIX-only empty-string off switch
	BOOST_CHECK(!bertini::records::AmbientRecordsPath().has_value());
	unsetenv("BERTINI_RECORDS_DIR");        // leave the environment as we found it
#else
	_putenv("BERTINI_RECORDS_DIR=");
#endif
}

BOOST_AUTO_TEST_CASE(annotate_convenience_appends_annotation_records)
{
	OutputDirectory out(FreshDir("annotate"));
	out.Annotate("abc123", 0, "projection", json::value(1.5));
	out.Annotate("abc123", 0, "projection", json::value(2.5));   // newest wins downstream
	auto const records = out.Scan();
	BOOST_REQUIRE_EQUAL(records.size(), 2u);
	BOOST_CHECK_EQUAL(std::string(records[1].at("kind").as_string()), "annotation");
	BOOST_CHECK_EQUAL(records[1].at("point").as_object().at("index").as_int64(), 0);
	BOOST_CHECK_EQUAL(records[1].at("value").as_double(), 2.5);
}

BOOST_AUTO_TEST_SUITE_END()
