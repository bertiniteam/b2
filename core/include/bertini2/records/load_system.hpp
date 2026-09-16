//This file is part of Bertini 2.
//
//load_system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//load_system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with load_system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file load_system.hpp

\brief Reloading an archived system from a records directory by its content digest.

The archive stores every system and homotopy a solve touched as a definition document
whose `encoding` member is the system's exact canonical encoding -- the digest preimage
(see SystemEncodingAsJson).  Loading is: fetch the document, decode the encoding
(System::FromCanonicalEncoding), check that the rebuilt system's digest is the digest
it was filed under, and intern it.  The digest check is the proof of fidelity: a loaded
system is the archived one, byte for byte in its identity, or the load fails.
*/

#pragma once

#include <memory>
#include <string>

#include "bertini2/records/output_directory.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {
namespace records {

/**
\brief The canonical encoding text inside an archived system definition document.

\param definition_json The document as PutDefinition stored it (the output of
       SystemEncodingAsJson): a JSON object with an `encoding` string member.
\return The encoding text, exactly as archived.

Throws std::runtime_error if the document is not JSON or has no `encoding` string.
*/
std::string EncodingOfSystemDefinition(std::string const& definition_json);

/**
\brief Load the system archived under `digest_hex` in `directory`.

\param directory The records directory holding the `systems/` definition.
\param digest_hex The system's content digest, as the run header and the definition
       filename carry it (64 lowercase hex characters).
\return The interned representative of the rebuilt system: if an equal system is already
        alive in this process, THAT one comes back.

Throws std::runtime_error if no definition has that id, if the encoding cannot be read
(a version this build does not read, other canonicalization settings, damage), or if the
rebuilt system's digest differs from `digest_hex` -- the document was altered, or the
encoding no longer means what it meant when it was written.
*/
std::shared_ptr<const System> LoadSystem(OutputDirectory const& directory, std::string const& digest_hex);

} // namespace records
} // namespace bertini
