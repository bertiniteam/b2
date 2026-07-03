//This file is part of Bertini 2.
//
//output_directory.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//output_directory.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with output_directory.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file output_directory.cpp

\brief The structured output directory implementation (ledgerrec/1, ADR-0045).
Boost.JSON is used header-only (src.hpp included here, exactly once in the library) so
no new link component is required on any platform.
*/

#include <boost/json/src.hpp>   // header-only Boost.JSON: this TU provides the impl

#include "bertini2/records/output_directory.hpp"

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace bertini {
namespace records {

namespace json = boost::json;

namespace {

	constexpr char kReadme[] =
R"(This directory is a structured output directory: the durable, self-contained
record of numerical algebraic geometry computations (polynomial-system solves by
homotopy continuation).  Written by bertini2 (record schema ledgerrec/1; full spec:
docs/records/ledgerrec-1.md in the bertini2 repository -- but this file suffices).
It needs no software to read, and you are free to delete it -- the only consequence
is recomputing.

LAYOUT
  results.json  the declared results, pretty-printed and SELF-COMPLETE: the final
                results first, then a "runs" section referring to everything used to
                construct them (system rendering, configs, seed) by definition id.
                Most readers start AND END here; it is one json.load away.
  INDEX.txt     one line per run: when, what was solved, how many paths.
  history/      the records: JSON, one object per line (JSONL), one file per writing
                session, named by date.  Read with eyes, grep, jq, or
                pandas.read_json(..., lines=True).
  definitions/  the things records refer to (polynomial systems, etc.), stored once
                each in a file named by the SHA-256 of its content -- so records can
                reference them exactly, and `sha256sum` verifies them.  Mostly
                human-readable input files.

RECORD FORMAT (schema ledgerrec/1) -- every history line is one JSON object:
  kind="run"        a solve: `ask` (what was requested: target system digest + config
                    + seed), `run` (this run's id), `when`, `num_paths`, and how start
                    points arise (recorded values, or a reference to an ancestor run).
  kind="track"      one continued path: `run`, `index`, `status`, `endpoint`
                    (coordinates as [real, imaginary] decimal-string pairs, full
                    precision), and `start` (its provenance: a start_label, or a
                    point_ref {run, index} into an ancestor run's endpoint).
  kind="result"     the declared DELIVERABLES: `name`, `points` [{run, index}, ...],
                    optional inline `value`.  Everything else in history/ is
                    scaffolding; results.json renders these -- "what were my solutions?".
  kind="annotation" metadata attached to a point: `point` {run, index}, `key`, `value`.
  kind="given"      externally supplied data: `source` (definition id) -- provenance
                    bottoms out honestly at the boundary of what was computed here.
Readers preserve records of kinds they do not recognize.  Chains of runs are walkable:
follow track records' `start` references backward until a start_label or a given --
that is the complete provenance of any point recorded here.
)";

	std::string TimeStamp(char const* fmt)
	{
		std::time_t now = std::time(nullptr);
		std::tm tm_buf{};
#ifdef _WIN32
		localtime_s(&tm_buf, &now);
#else
		localtime_r(&now, &tm_buf);
#endif
		char buffer[32];
		std::strftime(buffer, sizeof(buffer), fmt, &tm_buf);
		return buffer;
	}

	// A small pretty-printer (boost::json::serialize is compact-only): 1-space-indented,
	// newline-separated -- results.json is the file a human most interacts with.
	void PrettyPrint(std::ostream& out, json::value const& v, int depth)
	{
		std::string const pad(static_cast<std::size_t>(depth) + 1, ' ');
		std::string const pad_close(static_cast<std::size_t>(depth), ' ');
		switch (v.kind())
		{
			case json::kind::object:
			{
				auto const& obj = v.get_object();
				if (obj.empty()) { out << "{}"; return; }
				out << "{\n";
				bool first = true;
				for (auto const& kv : obj)
				{
					if (!first) out << ",\n";
					first = false;
					out << pad << json::serialize(json::value(kv.key())) << ": ";
					PrettyPrint(out, kv.value(), depth + 1);
				}
				out << "\n" << pad_close << "}";
				return;
			}
			case json::kind::array:
			{
				auto const& arr = v.get_array();
				if (arr.empty()) { out << "[]"; return; }
				// short leaf arrays (coordinate triples etc.) stay on one line
				bool leaf = true;
				for (auto const& e : arr)
					if (e.is_object() || e.is_array()) { leaf = false; break; }
				if (leaf && arr.size() <= 4) { out << json::serialize(v); return; }
				out << "[\n";
				bool first = true;
				for (auto const& e : arr)
				{
					if (!first) out << ",\n";
					first = false;
					out << pad;
					PrettyPrint(out, e, depth + 1);
				}
				out << "\n" << pad_close << "]";
				return;
			}
			default:
				out << json::serialize(v);
		}
	}

	std::string GetString(json::object const& obj, char const* key, std::string const& fallback = "")
	{
		auto const* v = obj.if_contains(key);
		if (v && v->is_string())
			return std::string(v->get_string());
		return fallback;
	}

} // unnamed namespace


OutputDirectory::OutputDirectory(std::filesystem::path root) : root_(std::move(root))
{
	std::filesystem::create_directories(root_ / "definitions");
	std::filesystem::create_directories(root_ / "history");
	auto const readme = root_ / "README.txt";
	if (!std::filesystem::exists(readme))
	{
		std::ofstream out(readme);
		out << kReadme;
	}
}


// ---- definitions ----

std::filesystem::path OutputDirectory::DefinitionPath(std::string const& id) const
{
	return root_ / "definitions" / id.substr(0, 2) / id.substr(2);
}

std::string OutputDirectory::PutDefinition(std::string const& content,
                                           std::optional<std::string> external_id)
{
	std::string const id = external_id ? *external_id : detail::Sha256(content).Hex();
	auto const path = DefinitionPath(id);
	if (std::filesystem::exists(path))
		return id;   // idempotent: content-addressed writes never conflict
	std::filesystem::create_directories(path.parent_path());
	auto const tmp = path.parent_path() / (path.filename().string() + ".tmp");
	{
		std::ofstream out(tmp, std::ios::binary);
		out.write(content.data(), static_cast<std::streamsize>(content.size()));
	}
	std::filesystem::rename(tmp, path);   // atomic on POSIX
	return id;
}

std::string OutputDirectory::GetDefinition(std::string const& id) const
{
	std::ifstream in(DefinitionPath(id), std::ios::binary);
	if (!in.good())
		throw std::runtime_error("OutputDirectory: no definition with id " + id);
	std::ostringstream contents;
	contents << in.rdbuf();
	return contents.str();
}

bool OutputDirectory::HasDefinition(std::string const& id) const
{
	return std::filesystem::exists(DefinitionPath(id));
}


// ---- history ----

void OutputDirectory::EnsureSessionFile()
{
	if (session_.is_open())
		return;
	auto const stamp = TimeStamp("%Y%m%d_%H%M%S");
	auto const base = root_ / "history";
#ifdef _WIN32
	auto const pid = static_cast<long>(_getpid());
#else
	auto const pid = static_cast<long>(getpid());
#endif
	std::string suffix;
	for (char c = 'b'; c <= 'z'; ++c)
	{
		auto path = base / (stamp + "-pid" + std::to_string(pid) + suffix + ".jsonl");
		// claim the name with an exclusive create ("x" mode): one writer per file, ever
		if (std::FILE* claimed = std::fopen(path.string().c_str(), "wx"))
		{
			std::fclose(claimed);
			session_ = std::ofstream(path, std::ios::app);
			session_path_ = path;
			return;
		}
		suffix = std::string(1, c);
	}
	throw std::runtime_error("OutputDirectory: could not claim a session history file");
}

void OutputDirectory::Append(json::object const& record)
{
	EnsureSessionFile();
	session_ << json::serialize(record) << "\n";
	session_.flush();
	auto const* kind = record.if_contains("kind");
	if (kind && kind->is_string() && kind->get_string() == "run")
		RefreshIndex();
}

std::vector<json::object> OutputDirectory::Scan() const
{
	std::vector<json::object> records;
	std::vector<std::filesystem::path> files;
	for (auto const& entry : std::filesystem::directory_iterator(root_ / "history"))
		if (entry.path().extension() == ".jsonl")
			files.push_back(entry.path());
	std::sort(files.begin(), files.end());

	for (auto const& path : files)
	{
		std::ifstream in(path);
		std::vector<std::string> lines;
		for (std::string line; std::getline(in, line); )
			lines.push_back(line);
		for (std::size_t ii = 0; ii < lines.size(); ++ii)
		{
			if (lines[ii].find_first_not_of(" \t\r") == std::string::npos)
				continue;
			boost::system::error_code ec;
			auto parsed = json::parse(lines[ii], ec);
			if (ec)
			{
				if (ii + 1 == lines.size())
					continue;   // torn tail: the write died mid-line; replay ignores it
				throw std::runtime_error("corrupt history file " + path.string()
					+ " at line " + std::to_string(ii + 1));
			}
			records.push_back(parsed.as_object());
		}
	}
	return records;
}


// ---- derived views ----

void OutputDirectory::RefreshIndex() const
{
	auto const records = Scan();
	std::map<std::string, std::pair<long, long>> track_counts;   // run -> (done, failed)
	for (auto const& r : records)
		if (GetString(r, "kind") == "track")
		{
			auto& counts = track_counts[GetString(r, "run")];
			++counts.first;
			if (GetString(r, "status") == "failed")
				++counts.second;
		}

	std::ostringstream out;
	out << "what has been solved here (newest last; details in history/):\n\n";
	for (auto const& r : records)
	{
		if (GetString(r, "kind") != "run")
			continue;
		auto const run_id = GetString(r, "run", "?");
		auto const counts = track_counts.count(run_id) ? track_counts[run_id] : std::make_pair(0L, 0L);
		std::string num_paths = "?";
		if (auto const* np = r.if_contains("num_paths"); np && np->is_int64())
			num_paths = std::to_string(np->get_int64());

		std::string description = "(unknown target)";
		auto const target_object = GetString(r, "target_object");
		if (!target_object.empty() && HasDefinition(target_object))
		{
			std::istringstream text(GetDefinition(target_object));
			std::string functions;
			for (std::string line; std::getline(text, line); )
			{
				auto const eq = line.find('=');
				if (eq == std::string::npos)
					continue;
				auto trimmed = line.substr(line.find_first_not_of(" \t"));
				while (!trimmed.empty() && (trimmed.back() == ';' || trimmed.back() == '\r' || trimmed.back() == ' '))
					trimmed.pop_back();
				if (!functions.empty())
					functions += "; ";
				functions += trimmed;
			}
			if (!functions.empty())
				description = functions.size() > 60 ? functions.substr(0, 57) + "..." : functions;
		}

		out << GetString(r, "when", "????-??-?? ??:??") << "  "
		    << counts.first << "/" << num_paths << " paths done";
		if (counts.second > 0)
			out << " (" << counts.second << " failed)";
		out << "  " << GetString(r, "op", "solve")
		    << "  " << description
		    << "   [run " << run_id << "]\n";
	}
	std::ofstream index(root_ / "INDEX.txt");
	index << out.str();
}

void OutputDirectory::RefreshResults() const
{
	auto const records = Scan();

	std::map<std::pair<std::string, std::int64_t>, json::object> tracks;
	std::map<std::pair<std::string, std::int64_t>, json::object> annotations;
	std::map<std::string, json::object> runs;
	std::map<std::string, json::object> declared;   // newest declaration of a name wins

	for (auto const& r : records)
	{
		auto const kind = GetString(r, "kind");
		if (kind == "run")
			runs[GetString(r, "run")] = r;
		else if (kind == "track")
		{
			if (auto const* idx = r.if_contains("index"); idx && idx->is_int64())
				tracks[{GetString(r, "run"), idx->get_int64()}] = r;
		}
		else if (kind == "annotation")
		{
			if (auto const* pt = r.if_contains("point"); pt && pt->is_object())
			{
				auto const& p = pt->get_object();
				if (auto const* idx = p.if_contains("index"); idx && idx->is_int64())
				{
					auto& bag = annotations[{GetString(p, "run"), idx->get_int64()}];
					bag[GetString(r, "key")] = r.contains("value") ? r.at("value") : json::value();
				}
			}
		}
		else if (kind == "result")
			declared[GetString(r, "name")] = r;
	}

	json::object results_section;
	json::object runs_section;   // only the runs the declared results reference

	for (auto const& [name, rec] : declared)
	{
		json::object entry;
		entry["declared"] = rec.contains("when") ? rec.at("when") : json::value();
		entry["description"] = GetString(rec, "description");
		if (rec.contains("value"))
			entry["value"] = rec.at("value");

		json::array points;
		if (auto const* refs = rec.if_contains("points"); refs && refs->is_array())
		{
			for (auto const& ref_value : refs->get_array())
			{
				if (!ref_value.is_object())
					continue;
				auto const& ref = ref_value.get_object();
				auto const run_id = GetString(ref, "run");
				std::int64_t index = -1;
				if (auto const* idx = ref.if_contains("index"); idx && idx->is_int64())
					index = idx->get_int64();

				// the runs section makes results.json SELF-COMPLETE: final results first,
				// then references to what constructed them (system, configs, seed)
				auto const run_it = runs.find(run_id);
				if (run_it != runs.end() && !runs_section.contains(run_id))
				{
					json::object summary;
					for (char const* key : {"when", "op", "ask", "target_object",
					                        "target_digest", "config_object", "num_paths"})
						if (run_it->second.contains(key))
							summary[key] = run_it->second.at(key);
					runs_section[run_id] = summary;
				}

				json::object point;
				point["provenance"] = ref;
				auto const track_it = tracks.find({run_id, index});
				if (track_it == tracks.end() || GetString(track_it->second, "status") != "success")
				{
					point["status"] = "missing";
					points.push_back(point);
					continue;
				}
				point["status"] = "success";

				// user variable names label coordinates only when the counts agree;
				// internal (homogenized) points have extra coordinates, and labeling
				// them with user names would be misleading (user-coordinate rendering
				// is a noted follow-up)
				std::vector<std::string> var_names;
				// the run header's "variables" field is authoritative: the solver's own
				// internal ordering, homogenizing variables included
				if (run_it != runs.end())
					if (auto const* vars = run_it->second.if_contains("variables");
					    vars && vars->is_array())
						for (auto const& v : vars->get_array())
							if (v.is_string())
								var_names.emplace_back(v.get_string());
				// older directories: fall back to parsing the classic rendering
				if (var_names.empty() && run_it != runs.end())
				{
					auto const target_object = GetString(run_it->second, "target_object");
					if (!target_object.empty() && HasDefinition(target_object))
					{
						std::istringstream text(GetDefinition(target_object));
						for (std::string line; std::getline(text, line); )
						{
							auto const key_pos = line.find("variable_group");
							auto const alt_pos = line.find("variable ");
							if (key_pos == std::string::npos && alt_pos == std::string::npos)
								continue;
							auto names = line.substr(line.find(' ') + 1);
							while (!names.empty() && (names.back() == ';' || names.back() == '\r'))
								names.pop_back();
							std::istringstream splitter(names);
							for (std::string v; std::getline(splitter, v, ','); )
							{
								auto const b = v.find_first_not_of(" \t");
								auto const e = v.find_last_not_of(" \t");
								if (b != std::string::npos)
									var_names.push_back(v.substr(b, e - b + 1));
							}
							break;
						}
					}
				}

				json::object coords;
				if (auto const* endpoint = track_it->second.if_contains("endpoint");
				    endpoint && endpoint->is_array())
				{
					if (var_names.size() != endpoint->get_array().size())
						var_names.clear();
					std::size_t k = 0;
					for (auto const& coordinate : endpoint->get_array())
					{
						std::string const var = k < var_names.size() ? var_names[k]
						                                             : ("coordinate_" + std::to_string(k));
						coords[var] = coordinate;
						++k;
					}
				}
				point["coordinates"] = coords;

				auto const note_it = annotations.find({run_id, index});
				point["annotations"] = (note_it != annotations.end()) ? note_it->second : json::object{};
				points.push_back(point);
			}
		}
		if (!points.empty())
			entry["points"] = points;
		results_section[name] = entry;
	}

	json::object machine;
	machine["results"] = results_section;
	machine["runs"] = runs_section;

	{
		std::ofstream out(root_ / "results.json");
		PrettyPrint(out, machine, 0);
		out << "\n";
	}
	// RESULTS.txt retired: it duplicated results.json, which is now pretty-printed and
	// self-complete -- one results file.  Remove a stale copy from older directories.
	std::error_code ignored;
	std::filesystem::remove(root_ / "RESULTS.txt", ignored);
}


std::string OutputDirectory::Describe() const
{
	long history_files = 0;
	for (auto const& entry : std::filesystem::directory_iterator(root_ / "history"))
		if (entry.path().extension() == ".jsonl")
			++history_files;
	long definitions = 0;
	for (auto const& entry : std::filesystem::recursive_directory_iterator(root_ / "definitions"))
		if (entry.is_regular_file())
			++definitions;
	std::ostringstream out;
	out << Scan().size() << " records in " << history_files << " history file(s), "
	    << definitions << " definition(s), at " << std::filesystem::absolute(root_).string();
	return out.str();
}

} // namespace records
} // namespace bertini
