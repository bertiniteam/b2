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

\brief The structured output directory implementation (b2rec/1, ADR-0045).
Boost.JSON is used header-only (src.hpp included here, exactly once in the library) so
no new link component is required on any platform.
*/

#include <boost/json/src.hpp>   // header-only Boost.JSON: this TU provides the impl

#include "bertini2/records/output_directory.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
homotopy continuation).  Written by bertini2 (record schema b2rec/1; full spec:
docs/records/b2rec-1.md in the bertini2 repository -- but this file suffices).
It needs no software to read, and you are free to delete it -- the only consequence
is recomputing.

LAYOUT -- three stores, separated by concern:
  history/      WHAT WAS ASKED, WHEN: run headers, declared results, annotations.
                JSON, one object per line (JSONL), one file per writing session,
                named by date.  Every line is small; the data it speaks of lives in
                the other two stores, referred to by id.  Read with eyes, grep, jq,
                or pandas.read_json(..., lines=True).
  results/      WHAT WAS COMPUTED: one append-only JSONL file per run, at
                    results/<2 hex>/<run id>.jsonl
                Line 1 is a self-description header {"kind":"results_header", "run",
                "ask"} -- a wandering file says what it is.  Then one line per
                completed path: {"kind":"path", "index", "status", endpoint
                coordinates in full precision, and its per-path metadata (verdict
                names, cycle number, timings)} -- the data and its facts travel
                together.  Paths append as they complete, so a partially-computed
                run reads honestly: what is here is done.
  definitions/  WHAT THINGS ARE: content-addressed inputs, filed as
                    definitions/<kind>/<2 hex>/<kind>[-<role>]-<digest>.<ext>
                e.g.  definitions/systems/03/system-03958a...7f.json
                      definitions/givens/34/given-cli_input-3468cd...9b.txt
                The filename carries the FULL digest (never concatenate) and the kind;
                the two-hex folder exists purely so no directory grows unbounded; the
                extension is honest (.json for JSON, .txt for text).  Kinds:
                  systems/  the exact polynomial systems -- targets AND the homotopies
                            actually tracked -- as JSON: {"schema", "digest",
                            "system", "encoding"}.  The "system" value has fields for
                            the parts (variable groups, path variable, functions,
                            patches).  The encoding is bertini2's canonical form
                            (b2sysenc; versioned, block structure preserved) and is
                            the digest PREIMAGE: the id equals the system's content
                            digest, and hashing the encoding
                            (`jq -r .encoding <file> | sha256sum`) reproduces it.
                  configs/  the solver settings that ran, as JSON (digest embedded).
                  givens/   externally supplied data: start points (JSON), CLI input
                            files (byte-exact copies of what you supplied --
                            `sha256sum` reproduces their id directly).
  INDEX.txt     derived, rebuildable: one line per run -- when, what was solved, how
                many paths done.

RECORD FORMAT (schema b2rec/1) -- every history line is one JSON object:
  kind="run"        a solve: `ask` (what was requested: target system digest +
                    homotopy digest + config digest + seed), `run` (this run's id),
                    `when`, `num_paths`, `results_file` (where its computed paths
                    live), and how start points arise (recorded values, or a
                    reference to an ancestor run).  `target_digest` is the solved
                    system's content digest; the definitions/systems/ file with that
                    id holds the system's parts and its canonical encoding
                    (dereferencing is the reader's job).
  kind="recall"     a re-ask answered from the store: `run`, `when`, `num_recalled`,
                    `num_computed`.  A point is computed exactly once, ever; "recalled"
                    is a property of a session, narrated here -- never marked on points.
  kind="result"     the declared DELIVERABLES: `name`, `points` [{run, index}, ...],
                    optional inline `value` -- "what were my solutions?".
  kind="annotation" metadata attached to a point: `point` {run, index}, `key`, `value`.
  kind="given"      externally supplied data: `source` (definition id) -- provenance
                    bottoms out honestly at the boundary of what was computed here.
And every results/ line: kind="results_header" (line 1), then kind="path" -- one
continued path: `index`, `status` (success / diverged / failed), `endpoint` and
`endpoint_user` (coordinates as [real, imaginary] decimal-string pairs, full
precision), per-path metadata, and `start` (its provenance: a start_label, or a
point_ref {run, index} into an ancestor run's endpoint).
Readers preserve records of kinds they do not recognize.  Chains of runs are walkable:
follow path records' `start` references backward until a start_label or a given --
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
	// One textual form for a double EVERYWHERE the records write one: the shortest
	// round-trip decimal (0.5, 1e-05), never Boost.JSON's uppercase-scientific (5E-1).
	// Keeps annotation values, config views, and every other double consistent.
	std::string ReadableDouble(double d)
	{
		if (!std::isfinite(d))
			return json::serialize(json::value(d));   // null, per JSON rules
		char buffer[32];
		auto const res = std::to_chars(buffer, buffer + sizeof(buffer), d);
		return std::string(buffer, res.ptr);
	}

	// json::serialize with ReadableDouble applied to every double leaf (compact form,
	// used for the one-line history records)
	std::string SerializeReadable(json::value const& v)
	{
		switch (v.kind())
		{
			case json::kind::double_:
				return ReadableDouble(v.get_double());
			case json::kind::object:
			{
				std::string out = "{";
				bool first = true;
				for (auto const& kv : v.get_object())
				{
					if (!first) out += ",";
					first = false;
					out += json::serialize(json::value(kv.key())) + ":" + SerializeReadable(kv.value());
				}
				return out + "}";
			}
			case json::kind::array:
			{
				std::string out = "[";
				bool first = true;
				for (auto const& e : v.get_array())
				{
					if (!first) out += ",";
					first = false;
					out += SerializeReadable(e);
				}
				return out + "]";
			}
			default:
				return json::serialize(v);
		}
	}

	void PrettyPrint(std::ostream& out, json::value const& v, int depth)
	{
		std::string const pad(static_cast<std::size_t>(depth) + 1, ' ');
		std::string const pad_close(static_cast<std::size_t>(depth), ' ');
		switch (v.kind())
		{
			case json::kind::double_:
				out << ReadableDouble(v.get_double());
				return;
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
				if (leaf && arr.size() <= 4) { out << SerializeReadable(v); return; }
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

	// Views (results.json, INDEX.txt) are REGENERATED files that concurrent writers may
	// refresh at once (e.g. two MPI ranks recording their own runs): write to a
	// process-unique temp then rename, so a reader never sees a truncated/empty view.
	// (Definitions already write this way; journals are one-writer-per-file.)
	void WriteViewAtomically(std::filesystem::path const& target, std::string const& content)
	{
		auto const tmp = target.parent_path() /
			(target.filename().string() + ".tmp." + std::to_string(
#ifdef _WIN32
				_getpid()
#else
				getpid()
#endif
			));
		{
			std::ofstream out(tmp, std::ios::binary);
			out.write(content.data(), static_cast<std::streamsize>(content.size()));
		}
		std::error_code ec;
		std::filesystem::rename(tmp, target, ec);   // atomic on POSIX; last writer wins
		if (ec)
			std::filesystem::remove(tmp, ec);
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
	std::filesystem::create_directories(root_ / "results");
	auto const readme = root_ / "README.txt";
	if (!std::filesystem::exists(readme))
	{
		std::ofstream out(readme);
		out << kReadme;
	}
}


std::shared_ptr<OutputDirectory> OutputDirectory::Shared(std::filesystem::path const& root)
{
	static std::mutex table_mutex;
	static std::map<std::filesystem::path, std::weak_ptr<OutputDirectory>> table;

	std::error_code ec;
	auto key = std::filesystem::weakly_canonical(root, ec);
	if (ec)
		key = root;

	std::lock_guard<std::mutex> lock(table_mutex);
	auto& slot = table[key];
	if (auto live = slot.lock())
		return live;
	auto made = std::make_shared<OutputDirectory>(root);
	slot = made;
	return made;
}

std::string TimeStampNow()
{
	return TimeStamp("%Y-%m-%d %H:%M");
}

std::optional<std::string> AmbientRecordsPath()
{
	char const* dir = std::getenv("BERTINI_RECORDS_DIR");
	if (dir == nullptr)
		return std::string("bertini_output");   // records on by default (ADR-0047)
	std::string const value(dir);
	if (value.empty() || value == "none")       // the explicit off switch ("" is POSIX-only)
		return std::nullopt;
	return value;
}


// ---- definitions ----

namespace {

	// "systems" -> "system": the filename carries the kind in the singular
	std::string KindSingular(std::string const& kind)
	{
		return (!kind.empty() && kind.back() == 's') ? kind.substr(0, kind.size() - 1) : kind;
	}

} // unnamed namespace

std::string SystemEncodingAsJson(std::string const& encoding_text, std::string const& digest_hex,
                                 json::object const& parts)
{
	// the schema token is the encoding's own first word, so the two never drift
	auto const schema_end = encoding_text.find_first_of(" \n");
	std::string const schema = encoding_text.substr(0, schema_end);
	// boost::json does the escaping (encodings may contain any variable name -- emoji
	// included); the parts object is pretty-printed so the file reads well, and the
	// encoding (one long string) goes last
	json::object doc;
	doc["schema"] = schema;
	doc["digest"] = digest_hex;
	doc["system"] = parts;
	doc["encoding"] = encoding_text;
	std::ostringstream out;
	PrettyPrint(out, doc, 0);
	out << "\n";
	return out.str();
}

namespace {

	// ids are 64 lowercase hex; kinds/labels are short lowercase words.  Anything else
	// is refused OUTRIGHT -- these strings become filesystem paths, and a hostile or
	// buggy caller must not be able to write outside definitions/ ("../../evil") or
	// smuggle separators into filenames.
	bool ValidDefinitionId(std::string const& id)
	{
		if (id.size() != 64)
			return false;
		for (char const c : id)
			if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
				return false;
		return true;
	}

	bool ValidPathWord(std::string const& word, bool allow_empty)
	{
		if (word.empty())
			return allow_empty;
		for (char const c : word)
			if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_'))
				return false;
		return true;
	}

} // unnamed namespace

std::filesystem::path OutputDirectory::DefinitionPath(std::string const& kind,
                                                      std::string const& id,
                                                      std::string const& label) const
{
	// kind folder for the browsing human; two-hex-char shard inside so no single
	// directory grows unbounded (a 100k-target sweep must not melt systems/).  The
	// filename repeats the FULL id -- recovering a definition's digest must never
	// require string concatenation -- plus an optional role label, so a listing says
	// what each file is even after it wanders away from its folder.  The id always
	// sits between the LAST '-' and the extension.
	return root_ / "definitions" / kind / id.substr(0, 2)
	       / (KindSingular(kind) + (label.empty() ? "" : "-" + label) + "-" + id);
	       // extension appended at write time
}

std::optional<std::filesystem::path> OutputDirectory::FindDefinition(std::string const& id) const
{
	// records carry bare ids; kind and extension are presentation.  Resolution: the
	// unique file under definitions/<kind>/<first 2 hex>/ named <kind>-<id>.<ext>.
	// Non-hex "ids" are refused: a hostile directory's records must not be able to
	// point resolution outside definitions/ (e.g. "../..").
	if (!ValidDefinitionId(id))
		return std::nullopt;
	auto const defs = root_ / "definitions";
	if (!std::filesystem::exists(defs))
		return std::nullopt;
	for (auto const& kind_dir : std::filesystem::directory_iterator(defs))
	{
		if (!kind_dir.is_directory())
			continue;
		auto const shard = kind_dir.path() / id.substr(0, 2);
		if (!std::filesystem::exists(shard))
			continue;
		for (auto const& entry : std::filesystem::directory_iterator(shard))
		{
			// filename shape: <kind>[-<label>]-<id>.<ext> -- the id is always
			// between the LAST '-' and the extension (hex never contains '-')
			auto const name = entry.path().filename().string();
			auto const dot = name.rfind('.');
			if (dot == std::string::npos)
				continue;
			auto const dash = name.rfind('-', dot);
			if (dash != std::string::npos && dot > dash
			    && name.compare(dash + 1, dot - dash - 1, id) == 0)
				return entry.path();
		}
	}
	return std::nullopt;
}

std::string OutputDirectory::PutDefinition(std::string const& content, std::string const& kind,
                                           std::optional<std::string> external_id,
                                           std::string const& label)
{
	if (!ValidPathWord(kind, false))
		throw std::invalid_argument("OutputDirectory: definition kind must be a lowercase "
		                            "word ([a-z0-9_]+), got '" + kind + "'");
	if (!ValidPathWord(label, true))
		throw std::invalid_argument("OutputDirectory: definition label must be a lowercase "
		                            "word ([a-z0-9_]*), got '" + label + "'");
	if (external_id && !ValidDefinitionId(*external_id))
		throw std::invalid_argument("OutputDirectory: external definition id must be 64 "
		                            "lowercase hex characters, got '" + *external_id + "'");
	std::string const id = external_id ? *external_id : detail::Sha256(content).Hex();
	if (auto const existing = FindDefinition(id))
		return id;   // idempotent: content-addressed writes never conflict
	char const* const ext = (!content.empty() && content.front() == '{') ? ".json" : ".txt";
	auto path = DefinitionPath(kind, id, label);
	path += ext;
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
	auto const path = FindDefinition(id);
	if (!path)
		throw std::runtime_error("OutputDirectory: no definition with id " + id);
	std::ifstream in(*path, std::ios::binary);
	std::ostringstream contents;
	contents << in.rdbuf();
	return contents.str();
}

bool OutputDirectory::HasDefinition(std::string const& id) const
{
	return FindDefinition(id).has_value();
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

void OutputDirectory::Annotate(std::string const& run_id, std::int64_t index,
                               std::string const& key, boost::json::value const& value)
{
	Append({{"kind", "annotation"},
	        {"point", json::object{{"run", run_id}, {"index", index}}},
	        {"key", key},
	        {"value", value}});
}

void OutputDirectory::Append(json::object const& record)
{
	std::lock_guard<std::mutex> lock(append_mutex_);
	EnsureSessionFile();
	session_ << SerializeReadable(record) << "\n";
	session_.flush();
	auto const* kind = record.if_contains("kind");
	if (kind && kind->is_string() && kind->get_string() == "run")
		RefreshIndex();
}

namespace {

	// The one JSONL reader, shared by history and results files: a torn FINAL line is
	// skipped (the crash-mid-append case), a torn interior line throws (real corruption).
	void ScanJsonlInto(std::filesystem::path const& path, std::vector<json::object>& records)
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
				throw std::runtime_error("corrupt records file " + path.string()
					+ " at line " + std::to_string(ii + 1));
			}
			records.push_back(parsed.as_object());
		}
	}

} // unnamed namespace

std::vector<json::object> OutputDirectory::Scan() const
{
	std::vector<json::object> records;
	std::vector<std::filesystem::path> files;
	for (auto const& entry : std::filesystem::directory_iterator(root_ / "history"))
		if (entry.path().extension() == ".jsonl")
			files.push_back(entry.path());
	std::sort(files.begin(), files.end());

	for (auto const& path : files)
		ScanJsonlInto(path, records);
	return records;
}


// ---- results (per-run payload files) ----

namespace {

	// run ids are hex (a prefix of the ask's SHA-256).  Anything else is refused: the
	// id becomes a filesystem path, and a hostile record must not write outside results/.
	bool ValidRunId(std::string const& id)
	{
		if (id.size() < 2 || id.size() > 64)
			return false;
		for (char const c : id)
			if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
				return false;
		return true;
	}

} // unnamed namespace

std::filesystem::path OutputDirectory::ResultsPath(std::string const& run_id) const
{
	// sharded like definitions/: two hex chars so no directory grows unbounded
	return root_ / "results" / run_id.substr(0, 2) / (run_id + ".jsonl");
}

void OutputDirectory::EnsureResultsFile(std::string const& run_id, json::object const& header)
{
	if (!ValidRunId(run_id))
		throw std::invalid_argument("OutputDirectory: invalid run id '" + run_id + "'");
	std::lock_guard<std::mutex> lock(append_mutex_);
	auto const path = ResultsPath(run_id);
	std::filesystem::create_directories(path.parent_path());
	// exclusive create: exactly one writer (across processes) writes the header line
	if (std::FILE* claimed = std::fopen(path.string().c_str(), "wx"))
	{
		std::fclose(claimed);
		std::ofstream out(path, std::ios::app);
		out << SerializeReadable(header) << "\n";
		out.flush();
	}
}

void OutputDirectory::AppendResult(std::string const& run_id, json::object const& record)
{
	if (!ValidRunId(run_id))
		throw std::invalid_argument("OutputDirectory: invalid run id '" + run_id + "'");
	std::lock_guard<std::mutex> lock(append_mutex_);
	auto& stream = results_streams_[run_id];
	if (!stream.is_open())
	{
		auto const path = ResultsPath(run_id);
		std::filesystem::create_directories(path.parent_path());
		stream.open(path, std::ios::app);
	}
	stream << SerializeReadable(record) << "\n";
	stream.flush();
}

std::vector<json::object> OutputDirectory::ResultsOf(std::string const& run_id) const
{
	std::vector<json::object> records;
	if (!ValidRunId(run_id))
		return records;
	auto const path = ResultsPath(run_id);
	if (std::filesystem::exists(path))
		ScanJsonlInto(path, records);
	return records;
}


// ---- derived views ----

void OutputDirectory::RefreshIndex() const
{
	auto const records = Scan();

	std::ostringstream out;
	out << "what has been solved here (newest last; paths in results/, details in history/):\n\n";
	for (auto const& r : records)
	{
		if (GetString(r, "kind") != "run")
			continue;
		auto const run_id = GetString(r, "run", "?");
		// paths done / failed, from the run's results file (the payload store).
		// Counted RAW -- lines and a byte pattern, never a JSON parse: this view
		// refreshes per solve, and a 300-million-path run must not be re-parsed to
		// render one INDEX line.  The pattern matches the compact serialization this
		// class itself writes; a torn tail undercounts by at most one (a view may).
		auto counts = std::make_pair(0L, 0L);
		if (ValidRunId(run_id))
			if (std::ifstream in(ResultsPath(run_id)); in)
				for (std::string line; std::getline(in, line); )
				{
					if (line.find("\"kind\":\"path\"") == std::string::npos)
						continue;
					++counts.first;
					if (line.find("\"status\":\"failed\"") != std::string::npos)
						++counts.second;
				}
		std::string num_paths = "?";
		if (auto const* np = r.if_contains("num_paths"); np && np->is_int64())
			num_paths = std::to_string(np->get_int64());

		std::string description = "(unknown target)";
		// the system definition's structured parts carry the functions; transitional
		// directories put a classic rendering in the header, and older ones stored the
		// rendering AS the definition -- read whichever this directory has
		std::string rendering = GetString(r, "target_rendering");
		auto target_object = GetString(r, "target_digest");
		if (target_object.empty())
			target_object = GetString(r, "target_object");   // older directories
		if (!target_object.empty() && HasDefinition(target_object))
		{
			auto const stored = GetDefinition(target_object);
			if (!stored.empty() && stored.front() == '{')
			{
				std::error_code parse_error;
				auto const doc = json::parse(stored, parse_error);
				if (!parse_error && doc.is_object())
					if (auto const* sys = doc.get_object().if_contains("system");
					    sys && sys->is_object())
						if (auto const* fns = sys->get_object().if_contains("functions");
						    fns && fns->is_array())
						{
							std::string joined;
							for (auto const& f : fns->get_array())
								if (f.is_string())
									joined += (joined.empty() ? "" : ";  ")
									          + std::string(f.get_string());
							if (!joined.empty())
								description = joined.size() > 60
								              ? joined.substr(0, 57) + "..." : joined;
						}
			}
			else if (rendering.empty())
				rendering = stored;   // oldest directories: the definition IS classic text
		}
		if (description == "(unknown target)" && !rendering.empty())
		{
			std::istringstream text(rendering);
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
	WriteViewAtomically(root_ / "INDEX.txt", out.str());
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
	long run_files = 0;
	if (std::filesystem::exists(root_ / "results"))
		for (auto const& entry : std::filesystem::recursive_directory_iterator(root_ / "results"))
			if (entry.is_regular_file() && entry.path().extension() == ".jsonl")
				++run_files;
	std::ostringstream out;
	out << Scan().size() << " records in " << history_files << " history file(s), "
	    << run_files << " run results file(s), "
	    << definitions << " definition(s), at " << std::filesystem::absolute(root_).string();
	return out.str();
}

} // namespace records
} // namespace bertini
