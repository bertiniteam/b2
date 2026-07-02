//This file is part of Bertini 2.
//
//derive.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//derive.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with derive.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file derive.hpp

\brief Seed-rooted, fully-specified random derivation (`b2rand/1`; ADR-0044; the
structured-output-directory arc, rung 2).

"Seed 42" is an identity: the same seed must produce the same draws forever — across
runs, platforms, standard libraries, and Boost versions.  std/boost random
DISTRIBUTIONS carry no such contract (their algorithms are implementation- or
release-dependent), so every identity-relevant draw in bertini2 routes through this
stream instead: SHA-256 in counter mode (the already-vendored, FIPS-pinned hash), with
draw algorithms specified exactly here and versioned as `b2rand/1`.  Any change to a
draw algorithm bumps the version — a new keyspace, never silent drift of what a seed
means.

Stream keying preserves the existing domain separation (setup / per-path / worker
streams; see random.hpp): a stream is keyed by SHA-256 over
`b2rand/1 || master seed || domain || index`, each as 8 big-endian bytes; block i of
the stream is SHA-256(key || i as 8 big-endian bytes).

Draw specifications (exact, forever):
- Uint64: the next 8 unconsumed stream bytes, big-endian.
- UnitDouble: (Uint64 >> 11) * 2^-53 — uniform on [0,1), exactly 53 bits.
- SymmetricDouble: 2*UnitDouble() - 1 — uniform on [-1,1] (both ops exact).
- Bits(n): the next ceil(n/8) stream bytes as a big-endian nonnegative integer, masked
  to n bits.
- IntSymmetric(B): rejection sampling — with n = bit-length of 2B, draw v = Bits(n)
  until v <= 2B; return v - B.  Uniform on [-B, B].
- UnitRealMp(digits): k = ceil(digits * log2(10)) + 1 bits; value = ldexp(M, -k) where
  M = Bits(k), materialized at `digits` decimal digits of precision.  The mpz→mpfr
  conversion and the power-of-two scaling are exact/correctly-rounded, so the value is
  bit-identical everywhere.
*/

#pragma once

#include <cstdint>

#include "bertini2/detail/sha256.hpp"
#include "bertini2/mpfr_extensions.hpp"   // mpz_int / real_mp (NOT num_traits: that would cycle through random.hpp)
#include "bertini2/records/draw_functions.hpp"

namespace bertini {
namespace records {

/// \brief The version tag of the draw-derivation scheme; bump on any algorithm change.
constexpr char RandomDerivationVersion[] = "b2rand/1";

/**
\brief A deterministic, fully-specified random stream (SHA-256 counter mode).

One per thread in practice (see ThreadDrawStream); reseeded by the same
(master, domain, index) tuples the legacy engine seeding uses, so setup / per-path /
worker domain separation is preserved.
*/
class DrawStream
{
public:
	/// \brief Construct unseeded; the first draw self-seeds from entropy (unseeded runs
	/// stay random).  Seeded runs call Reseed before any draw.
	DrawStream() = default;

	/// \brief Rekey the stream from (master seed, domain tag, stream index) under
	/// `b2rand/1` and restart its counter.
	void Reseed(std::uint64_t master, std::uint64_t domain, std::uint64_t index);

	/// \brief The next 8 stream bytes as a big-endian unsigned 64-bit integer.
	std::uint64_t Uint64();

	/// \brief Uniform on [0,1), exactly 53 bits: (Uint64 >> 11) * 2^-53.
	double UnitDouble();

	/// \brief Uniform on [-1,1]: 2*UnitDouble() - 1 (both operations exact).
	double SymmetricDouble();

	/// \brief The next n bits as a nonnegative integer (big-endian bytes, masked).
	mpz_int Bits(unsigned num_bits);

	/// \brief Uniform integer on [-bound, bound], by rejection over Bits(bitlen(2*bound)).
	mpz_int IntSymmetric(mpz_int const& bound);

	/// \brief Uniform on [0,1) at `digits` decimal digits: ldexp(Bits(k), -k) with
	/// k = ceil(digits*log2(10)) + 1, materialized at `digits` digits of precision.
	real_mp UnitRealMp(unsigned digits);

private:
	// Refill block_ with SHA-256(key || counter) and advance the counter.
	void NextBlock();
	// The next unconsumed stream byte (self-seeding from entropy if never seeded).
	unsigned char NextByte();

	detail::Digest256 key_;      ///< The stream key (hash of version, master, domain, index).
	std::uint64_t counter_ = 0;  ///< The block counter.
	unsigned char block_[32];    ///< The current block of stream bytes.
	unsigned offset_ = 32;       ///< Consumption offset into block_ (32 = exhausted).
	bool seeded_ = false;        ///< Whether Reseed (or entropy self-seeding) has run.
};

/**
\brief This thread's pinned draw stream — the source of every identity-relevant random
draw in bertini2.

Reseeded alongside the legacy engine by SetGlobalSeed / ReseedThisThread (random.hpp),
preserving the setup / per-path / worker stream separation.
*/
DrawStream& ThreadDrawStream();

} // namespace records
} // namespace bertini
