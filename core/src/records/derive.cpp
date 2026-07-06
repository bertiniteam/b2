//This file is part of Bertini 2.
//
//derive.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//derive.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with derive.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file derive.cpp

\brief The pinned draw stream (ADR-0044).  Every byte and every algorithm here is part
of the b2rand/1 forever-contract: changes bump the version, never silently alter what a
seed means.
*/

#include "bertini2/records/derive.hpp"

#include <random>
#include <string>

namespace bertini {
namespace records {

namespace {

	void AppendBigEndian64(std::string& s, std::uint64_t v)
	{
		for (int shift = 56; shift >= 0; shift -= 8)
			s.push_back(static_cast<char>((v >> shift) & 0xff));
	}

} // unnamed namespace


void DrawStream::Reseed(std::uint64_t master, std::uint64_t domain, std::uint64_t index)
{
	std::string key_material(RandomDerivationVersion);
	AppendBigEndian64(key_material, master);
	AppendBigEndian64(key_material, domain);
	AppendBigEndian64(key_material, index);
	key_ = detail::Sha256(key_material);
	counter_ = 0;
	offset_ = 32;
	seeded_ = true;
}

void DrawStream::NextBlock()
{
	std::string block_material(reinterpret_cast<char const*>(key_.bytes.data()), key_.bytes.size());
	AppendBigEndian64(block_material, counter_);
	++counter_;
	auto const digest = detail::Sha256(block_material);
	for (unsigned ii = 0; ii < 32; ++ii)
		block_[ii] = digest.bytes[ii];
	offset_ = 0;
}

unsigned char DrawStream::NextByte()
{
	if (!seeded_)
	{
		// unseeded runs stay random: self-seed from entropy (a seeded run has already
		// called Reseed via SetGlobalSeed / ReseedThisThread before any draw)
		std::random_device rd;
		Reseed((static_cast<std::uint64_t>(rd()) << 32) ^ rd(), 0, 0);
	}
	if (offset_ >= 32)
		NextBlock();
	return block_[offset_++];
}

std::uint64_t DrawStream::Uint64()
{
	std::uint64_t v = 0;
	for (unsigned ii = 0; ii < 8; ++ii)
		v = (v << 8) | NextByte();
	return v;
}

double DrawStream::UnitDouble()
{
	// exactly 53 bits: uniform on [0,1), every value representable
	return static_cast<double>(Uint64() >> 11) * 0x1.0p-53;
}

double DrawStream::SymmetricDouble()
{
	return 2.0 * UnitDouble() - 1.0;   // both operations exact
}

mpz_int DrawStream::Bits(unsigned num_bits)
{
	unsigned const num_bytes = (num_bits + 7) / 8;
	mpz_int v = 0;
	for (unsigned ii = 0; ii < num_bytes; ++ii)
		v = (v << 8) | NextByte();
	unsigned const excess = num_bytes * 8 - num_bits;
	return v >> excess;
}

mpz_int DrawStream::IntSymmetric(mpz_int const& bound)
{
	mpz_int const span = 2 * bound;
	unsigned num_bits = 0;
	for (mpz_int probe = span; probe > 0; probe >>= 1)
		++num_bits;
	while (true)
	{
		mpz_int const v = Bits(num_bits);
		if (v <= span)
			return v - bound;
	}
}

real_mp DrawStream::UnitRealMp(unsigned digits)
{
	// k = ceil(digits * log2(10)) + 1, computed in integer arithmetic:
	// log2(10) = 3.32192809...; 3.32193 > log2(10), so ceil via (digits*332193 + 99999)/100000
	unsigned const k = static_cast<unsigned>((static_cast<std::uint64_t>(digits) * 332193u + 99999u) / 100000u) + 1u;
	mpz_int const mantissa = Bits(k);
	// materialize AT the target precision (constructing at default precision first would
	// round there and lose bits); the mpz->mpfr conversion is correctly rounded, and the
	// power-of-two scaling is exact -- bit-identical on every platform
	real_mp value(mantissa, digits);
	using boost::multiprecision::ldexp;
	return ldexp(value, -static_cast<int>(k));
}


DrawStream& ThreadDrawStream()
{
	thread_local DrawStream stream;
	return stream;
}


// ---- the dependency-free draw functions (declared in draw_functions.hpp) ----

double DrawUnitDouble()
{
	return ThreadDrawStream().UnitDouble();
}

double DrawSymmetricDouble()
{
	return ThreadDrawStream().SymmetricDouble();
}

} // namespace records
} // namespace bertini
