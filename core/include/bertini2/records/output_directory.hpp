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
computations (`ledgerrec/1`; ADR-0045; the arc's rung 3).

Plain files are the source of truth, interactable without special software:
`history/` holds append-only JSONL records (one date-named file per writing session,
one writer per file, torn-final-line tolerant), `definitions/` holds content-addressed
definitions (atomic + idempotent writes; no locks exist or are needed), and README /
INDEX / RESULTS / results.json are derived, rebuildable views.  The format contract
lives in docs/records/ledgerrec-1.md and travels inside every directory as its
README.txt.  Cross-implementation compatibility with the Python pilot
(prototypes/ledger_v0) is tested.
*/

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include <boost/json.hpp>

#include "bertini2/detail/sha256.hpp"

namespace bertini {
namespace records {

/// \brief The record-format version tag written on run headers; see
/// docs/records/ledgerrec-1.md.  Bump on any change to record shapes or conventions.
constexpr char RecordSchemaVersion[] = "ledgerrec/1";

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

	/// \brief The directory's root path.
	std::filesystem::path const& Root() const { return root_; }

	// ---- definitions (content-addressed) ----

	/**
	\brief Store a definition content-addressed; returns its id.

	\param content The definition's bytes (usually human-readable text).
	\param external_id If given (e.g. a System's ContentDigest hex, whose preimage is
	       the canonical encoding rather than these bytes), store under that id;
	       otherwise the id is the SHA-256 of the bytes (self-verifying).
	\return The definition id (64 lowercase hex characters).

	Atomic (write-temp + rename) and idempotent: equal content lands at an equal path,
	so concurrent writers race benignly.
	*/
	std::string PutDefinition(std::string const& content,
	                          std::optional<std::string> external_id = std::nullopt);

	/// \brief Read a definition's bytes by id.  Throws if absent.
	std::string GetDefinition(std::string const& id) const;

	/// \brief Whether a definition with this id is present.
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
	\brief Every record from every history file, in filename-then-line order.

	A torn FINAL line of a file is skipped (the crash-mid-append case); a torn interior
	line throws (real corruption).  Records with unrecognized kinds are preserved.
	*/
	std::vector<boost::json::object> Scan() const;

	// ---- derived views (rebuildable at will) ----

	/// \brief (Re)write INDEX.txt: one line per run — when, paths done, op, target
	/// description, run id.  Derived from the records; never crashes a solve (defensive
	/// against malformed records).
	void RefreshIndex() const;

	/// \brief (Re)write results.json -- pretty-printed and self-complete -- from the
	/// declared `result` records: named results with coordinates keyed by variable
	/// name, annotations, provenance refs, and inline saved values.
	void RefreshResults() const;

	/// \brief One human line: how much is here (record/file/definition counts + path).
	std::string Describe() const;

private:
	std::filesystem::path DefinitionPath(std::string const& id) const;
	void EnsureSessionFile();

	std::filesystem::path root_;        ///< The directory root.
	std::ofstream session_;             ///< This session's history file (open after first append).
	std::filesystem::path session_path_; ///< Path of the session history file (empty until claimed).
};

} // namespace records
} // namespace bertini
