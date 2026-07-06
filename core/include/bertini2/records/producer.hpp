//This file is part of Bertini 2.
//
//producer.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//producer.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with producer.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file producer.hpp

\brief Which software produced a record: name, version, and source commit, for the
`producer` field of run headers.

An archive read in ten years should say which binary wrote it.  The producer is
DESCRIPTIVE only -- it is never part of the ask identity, so a newer build answering
the same ask still recalls the old answer rather than recomputing.
*/

#pragma once

#include <string>

#include <boost/json.hpp>

namespace bertini {
namespace records {

/// \brief This build's package version (e.g. "3.0.0.dev7", from the VERSION file).
std::string ProducerVersion();

/// \brief The git commit this build came from ("<sha>", "<sha>-dirty", or "unknown"
/// for builds outside a git checkout, e.g. from a source tarball).
std::string ProducerCommit();

/**
\brief The producer object recorded in every run header.

\return `{"name": "bertini2", "version": <version>, "commit": <sha or "unknown">}`.
*/
boost::json::object ProducerInfo();

} // namespace records
} // namespace bertini
