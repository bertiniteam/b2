//This file is part of Bertini 2.
//
//sha256.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//sha256.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with sha256.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file sha256.hpp

\brief Self-contained SHA-256 (FIPS 180-4) and the Digest256 value type.

This exists to give Systems a PERSISTENT content identity (ADR-0042): a digest that is
bit-identical across runs, compilers, standard libraries, and Boost versions.  That rules out
std::hash (implementation-defined), typeid().hash_code() (per-process), and
boost::uuids::detail::sha1 (detail:: namespace = no stability contract).  So we carry our own
~150-line clean-room FIPS 180-4 implementation, pinned by FIPS test vectors in
sha256_test.cpp.  It hashes canonical-encoding text (small, cold-path) -- throughput is
irrelevant; stability is everything.
*/

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <string_view>

namespace bertini {
namespace detail {

/**
\brief A 256-bit content digest -- the output of Sha256, used as a persistent identity key.

Value type with equality and ordering, so it can key std::map (the System intern table) and
be compared across sessions via its hex rendering.
*/
struct Digest256
{
	/// The raw 32 digest bytes, big-endian word order per FIPS 180-4.
	std::array<std::uint8_t, 32> bytes{};

	/// \brief The digest as 64 lowercase hex characters (the cross-session interchange form).
	std::string Hex() const;

	/// \brief Bytewise equality.
	bool operator==(Digest256 const& other) const { return bytes == other.bytes; }
	/// \brief Bytewise inequality.
	bool operator!=(Digest256 const& other) const { return bytes != other.bytes; }
	/// \brief Lexicographic byte order, so Digest256 can key ordered containers.
	bool operator<(Digest256 const& other) const { return bytes < other.bytes; }
};

/**
\brief SHA-256 (FIPS 180-4) of a byte string.

\param data The bytes to hash (canonical-encoding text, in practice).
\return The 256-bit digest.
*/
Digest256 Sha256(std::string_view data);

} // namespace detail
} // namespace bertini
