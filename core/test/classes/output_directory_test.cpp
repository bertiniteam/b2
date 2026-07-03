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
\file Tests for the C++ structured output directory (ledgerrec/1, ADR-0045): ports of
the Python pilot's ledger tests, plus the cross-implementation bridge — this suite
WRITES an example directory (into the ctest working directory) that the Python
prototype's cross-impl test reads, and READS one the prototype writes when present.
*/

#include <filesystem>
#include <fstream>

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
	BOOST_CHECK(text.find("ledgerrec/1") != std::string::npos);
	BOOST_CHECK(text.find("free to delete") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(definitions_are_idempotent_and_self_verifying)
{
	OutputDirectory out(FreshDir("defs"));
	auto const id1 = out.PutDefinition("hello records");
	auto const id2 = out.PutDefinition("hello records");
	BOOST_CHECK_EQUAL(id1, id2);
	BOOST_CHECK_EQUAL(out.GetDefinition(id1), "hello records");
	BOOST_CHECK_EQUAL(id1, bertini::detail::Sha256("hello records").Hex());

	auto const ext = out.PutDefinition("a rendering", std::string(64, 'a'));
	BOOST_CHECK_EQUAL(ext, std::string(64, 'a'));
	BOOST_CHECK(out.HasDefinition(ext));
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

BOOST_AUTO_TEST_CASE(index_and_results_render)
{
	auto dir = FreshDir("views");
	OutputDirectory out(dir);
	auto const target_id = out.PutDefinition(
		"function f;\nvariable_group x, y;\nf = x^2+4*y^2-4;\n");

	out.Append({{"kind", "run"}, {"schema", "ledgerrec/1"}, {"when", "2026-07-03 10:00"},
	            {"run", "abc123"}, {"op", "solve"}, {"target_object", target_id},
	            {"num_paths", 1}});
	out.Append({{"kind", "track"}, {"run", "abc123"}, {"index", 0}, {"status", "success"},
	            {"endpoint", json::array{json::array{"1.5", "0"}, json::array{"-0.66", "0"}}}});
	out.Append({{"kind", "annotation"}, {"point", json::object{{"run", "abc123"}, {"index", 0}}},
	            {"key", "projection"}, {"value", 1.5}});
	out.Append({{"kind", "result"}, {"name", "my solutions"}, {"description", "demo"},
	            {"when", "2026-07-03 10:01"},
	            {"points", json::array{json::object{{"run", "abc123"}, {"index", 0}}}}});
	out.RefreshResults();

	std::ifstream index(dir / "INDEX.txt");
	std::string index_text((std::istreambuf_iterator<char>(index)), std::istreambuf_iterator<char>());
	BOOST_CHECK(index_text.find("x^2+4*y^2-4") != std::string::npos);
	BOOST_CHECK(index_text.find("abc123") != std::string::npos);

	auto const results = json::parse([&]{
		std::ifstream in(dir / "results.json");
		return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	}()).as_object();
	auto const& mine = results.at("results").as_object().at("my solutions").as_object();
	auto const& point = mine.at("points").as_array()[0].as_object();
	BOOST_CHECK_EQUAL(std::string(point.at("status").as_string()), "success");
	BOOST_CHECK(point.at("coordinates").as_object().contains("x"));
	BOOST_CHECK(point.at("coordinates").as_object().contains("y"));
	BOOST_CHECK_EQUAL(point.at("annotations").as_object().at("projection").as_double(), 1.5);

	// self-completeness: the runs section refers back to what constructed the result
	auto const& run = results.at("runs").as_object().at("abc123").as_object();
	BOOST_CHECK_EQUAL(std::string(run.at("target_object").as_string()), target_id);

	// pretty-printed (the file a human most interacts with), and RESULTS.txt is gone
	{
		std::ifstream in(dir / "results.json");
		std::string raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		BOOST_CHECK(raw.find('\n') != std::string::npos);
	}
	BOOST_CHECK(!fs::exists(dir / "RESULTS.txt"));
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
