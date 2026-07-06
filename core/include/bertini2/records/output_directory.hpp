//This file is part of Bertini 2.
//
//output_directory.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//output_directory.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with output_directory.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file output_directory.hpp

\brief The structured output directory: durable, self-documenting records of
computations (`b2rec/1`; ADR-0045; the arc's rung 3).

Plain files are the source of truth, interactable without special software.  The three
truth stores separate concerns: `definitions/` holds the INPUTS by content address
(atomic + idempotent writes; no locks exist or are needed); `results/` holds the
OUTPUTS -- one append-only JSONL file per run, a self-description header line then one
line per completed path, everything recall needs; `history/` is the narrative of WHAT
WAS ASKED, WHEN -- run headers, result declarations, annotations -- every line small,
referring into the other two stores by id (one date-named file per writing session,
one writer per file, torn-final-line tolerant).  README / INDEX are derived,
rebuildable text views.  The format contract lives in docs/records/b2rec-1.md and
travels inside every directory as its README.txt.
*/

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include <boost/json.hpp>

#include "bertini2/detail/sha256.hpp"

namespace bertini {
namespace records {

/// \brief The record-format version tag written on run headers; see
/// docs/records/b2rec-1.md.  Bump on any change to record shapes or conventions.
constexpr char RecordSchemaVersion[] = "b2rec/1";

/**
\brief One structured output directory: definitions/ + history/ + derived views.

Construction ensures the directory skeleton and its self-documenting README.txt exist.
One OutputDirectory instance = one writing session = at most one history file (claimed
exclusively on first append).
*/
class OutputDirectory
{
public:
	/// \brief Open (creating if needed) the output directory at `root`, writing the
	/// self-documenting README.txt on first creation.
	explicit OutputDirectory(std::filesystem::path root);

	/**
	\brief The process-shared OutputDirectory at `root`: one instance -- hence one
	session history file -- per directory per process, however many solvers attach.

	Each instance claims its own session history file on first append, so anything
	constructing many writers in a loop (a parameter sweep attaching the ambient
	records to every solver, say) must share one instance: the per-second claim
	namespace is finite, and the history should not be shredded across a file per
	solve.  Keyed by the weakly-canonical path; instances are held weakly, so a
	directory nobody references anymore is released.

	\param root The directory root (need not exist yet).
	\return The shared instance for this path in this process.
	*/
	static std::shared_ptr<OutputDirectory> Shared(std::filesystem::path const& root);

	/// \brief The directory's root path.
	std::filesystem::path const& Root() const { return root_; }

	// ---- definitions (content-addressed, grouped by kind) ----

	/**
	\brief Store a definition content-addressed; returns its id.

	\param content The definition's bytes (usually human-readable text).
	\param kind The human-navigable subfolder the definition lives in: "systems",
	       "configs", or "givens" (an open vocabulary -- any lowercase name works;
	       these three are the ones bertini writes).  Inside each kind, files shard
	       by the id's first two hex characters so no directory grows unbounded.
	\param external_id If given, store under that id; otherwise the id is the SHA-256
	       of the bytes (self-verifying: `sha256sum` of the file reproduces its name).
	\param label Optional role label woven into the filename (e.g. "cli_input" gives
	       `given-cli_input-<digest>.txt`), so a listing says what each file IS
	       without opening it.  Presentation only, never identity.
	\return The definition id (64 lowercase hex characters).

	Atomic (write-temp + rename) and idempotent: equal content lands at an equal path,
	so concurrent writers race benignly.  Ids are resolved WITHOUT the kind (records
	reference bare ids); the kind is presentation, not identity.
	*/
	std::string PutDefinition(std::string const& content, std::string const& kind,
	                          std::optional<std::string> external_id = std::nullopt,
	                          std::string const& label = {});

	/// \brief Read a definition's bytes by id (searched across all kinds).  Throws if absent.
	std::string GetDefinition(std::string const& id) const;

	/// \brief Whether a definition with this id is present (searched across all kinds).
	bool HasDefinition(std::string const& id) const;

	// ---- history (append-only JSONL) ----

	/**
	\brief Append one record to this session's history file.

	The file is claimed on first append with an exclusive-create, date-named
	`YYYYMMDD_HHMMSS-pid<pid>[suffix].jsonl` (one writer per file, ever); each record is
	one compact JSON line, flushed per append so a kill loses at most the line being
	written (which Scan tolerates).  Appending a run header refreshes INDEX.txt.
	*/
	void Append(boost::json::object const& record);

	/**
	\brief Attach metadata to a recorded point: appends an annotation record
	`{"kind":"annotation", "point":{"run","index"}, "key", "value"}`.

	Annotations are the audit trail's margin notes -- projection values, "this is the
	one I meant", classification flags.  Readers (e.g. `bertini.load`) merge them
	beside the point they describe.  Newest wins per (point, key).

	\param run_id The run the point belongs to.
	\param index The point's path index within that run.
	\param key The annotation's name.
	\param value Any JSON value.
	*/
	void Annotate(std::string const& run_id, std::int64_t index,
	              std::string const& key, boost::json::value const& value);

	/**
	\brief Every record from every history file, in filename-then-line order.

	A torn FINAL line of a file is skipped (the crash-mid-append case); a torn interior
	line throws (real corruption).  Records with unrecognized kinds are preserved.
	*/
	std::vector<boost::json::object> Scan() const;

	// ---- results (per-run append-only JSONL payload files; truth) ----

	/**
	\brief Ensure the run's results file exists, with its self-description header as
	line 1: `{"kind":"results_header","schema",...,"run","ask"}`.

	The file lives at `results/<2 hex>/<run id>.jsonl` (sharded like definitions/).
	Creation is exclusive-create, so concurrent writers race benignly: exactly one
	writes the header.  Idempotent -- an existing file (a resumed run) is untouched.
	Written BEFORE the history run header that refers to it, so a reference never
	dangles.

	\param run_id The run id (lowercase hex; anything else is refused -- it becomes a path).
	\param header The self-description object written as line 1.
	*/
	void EnsureResultsFile(std::string const& run_id, boost::json::object const& header);

	/**
	\brief Append one payload record to the run's results file (flushed per record, so
	a kill loses at most the line in flight).

	Tracked paths append `{"kind":"path","run","index","status",...}` lines carrying
	the endpoint AND its per-path metadata -- the treasure travels with its facts.  The
	kind vocabulary is open: future operations may append other payload kinds.

	\param run_id The run whose file receives the record.
	\param record The payload record.
	*/
	void AppendResult(std::string const& run_id, boost::json::object const& record);

	/**
	\brief Every record in the run's results file, header line included, in order.

	Same tolerance rules as Scan: a torn final line is skipped, a torn interior line
	throws.  An absent file (nothing recorded for this run) is an empty vector.

	\param run_id The run to read.
	\return The records, oldest first; entry 0 is normally the results_header.
	*/
	std::vector<boost::json::object> ResultsOf(std::string const& run_id) const;

	// ---- derived views (rebuildable at will) ----

	/// \brief (Re)write INDEX.txt: one line per run — when, paths done, op, target
	/// description, run id.  Derived from the records; never crashes a solve (defensive
	/// against malformed records).
	void RefreshIndex() const;

	/// \brief One human line: how much is here (record/file/definition counts + path).
	std::string Describe() const;

private:
	std::filesystem::path DefinitionPath(std::string const& kind, std::string const& id,
	                                     std::string const& label) const;
	std::optional<std::filesystem::path> FindDefinition(std::string const& id) const;
	std::filesystem::path ResultsPath(std::string const& run_id) const;
	void EnsureSessionFile();

	std::filesystem::path root_;        ///< The directory root.
	std::ofstream session_;             ///< This session's history file (open after first append).
	std::filesystem::path session_path_; ///< Path of the session history file (empty until claimed).
	std::mutex append_mutex_;           ///< Serializes appends: a Shared() instance may be written from several threads.
	std::map<std::string, std::ofstream> results_streams_;  ///< Open per-run results files (keyed by run id; bounded by concurrent runs).
};

/**
\brief The archived form of a system definition: a machine-parseable JSON document
`{"schema", "digest", "system", "encoding"}` -- one `json.load` away, like the
config definitions.

The `system` value is the structured parts view (variable groups, path variable,
functions, patches -- see `io::SystemPartsJson`); the `encoding` value is the EXACT
canonical encoding text (`b2sysenc/<n>`, the digest preimage, block structure
preserved).  Classic (Bertini 1) syntax appears nowhere: it is an input/compatibility
format, not an output format.  Verification of a wandering file: extract the encoding
and hash it, e.g. `jq -r .encoding <file> | sha256sum` reproduces `digest`.

\param encoding_text The system's canonical encoding (`System::CanonicalEncodingText()`).
\param digest_hex The system's content digest (64 lowercase hex characters).
\param parts The structured parts view of the same system (presentation only).
\return The pretty-printed JSON document.
*/
std::string SystemEncodingAsJson(std::string const& encoding_text,
                                 std::string const& digest_hex,
                                 boost::json::object const& parts);

/// \brief The record timestamp, now: `YYYY-MM-DD HH:MM` local time -- the one clock
/// format every record kind uses (run headers, declarations, recall events).
std::string TimeStampNow();

/**
\brief Resolve the ambient records directory from the environment -- the ONE place the
`BERTINI_RECORDS_DIR` semantics live (the CLI and every solver's ambient attach use it).

Records are ON BY DEFAULT for every face of bertini2 (ADR-0047): when the variable is
unset, the resolution is `bertini_output`.  A non-empty value chooses the directory.
The OFF switch is the value `none` -- or the empty string, which only exists on POSIX
(Windows cannot represent an empty environment value: assigning one deletes the
variable, which would silently flip "off" back to the default).

\return The directory to record to, or nullopt when records are explicitly off.
*/
std::optional<std::string> AmbientRecordsPath();

} // namespace records
} // namespace bertini
