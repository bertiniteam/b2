//This file is part of Bertini 2.
//
//sha256.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//sha256.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with sha256.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file sha256.cpp

\brief Clean-room SHA-256 per FIPS 180-4 (single-shot, byte-oriented messages).
*/

#include "bertini2/detail/sha256.hpp"

#include <cstring>

namespace bertini {
namespace detail {

namespace {

// The FIPS 180-4 SHA-256 round constants: first 32 bits of the fractional parts of the cube
// roots of the first 64 primes.
constexpr std::uint32_t K[64] = {
	0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
	0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
	0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
	0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
	0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
	0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
	0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
	0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

inline std::uint32_t Rotr(std::uint32_t x, unsigned n) { return (x >> n) | (x << (32 - n)); }

// One 64-byte block through the compression function, updating the running hash state h.
void Compress(std::uint32_t h[8], std::uint8_t const block[64])
{
	std::uint32_t w[64];
	for (unsigned t = 0; t < 16; ++t)
		w[t] = (std::uint32_t(block[4*t]) << 24) | (std::uint32_t(block[4*t+1]) << 16)
		     | (std::uint32_t(block[4*t+2]) << 8) |  std::uint32_t(block[4*t+3]);
	for (unsigned t = 16; t < 64; ++t)
	{
		std::uint32_t const s0 = Rotr(w[t-15], 7) ^ Rotr(w[t-15], 18) ^ (w[t-15] >> 3);
		std::uint32_t const s1 = Rotr(w[t-2], 17) ^ Rotr(w[t-2], 19) ^ (w[t-2] >> 10);
		w[t] = w[t-16] + s0 + w[t-7] + s1;
	}

	std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
	for (unsigned t = 0; t < 64; ++t)
	{
		std::uint32_t const S1 = Rotr(e, 6) ^ Rotr(e, 11) ^ Rotr(e, 25);
		std::uint32_t const ch = (e & f) ^ (~e & g);
		std::uint32_t const temp1 = hh + S1 + ch + K[t] + w[t];
		std::uint32_t const S0 = Rotr(a, 2) ^ Rotr(a, 13) ^ Rotr(a, 22);
		std::uint32_t const maj = (a & b) ^ (a & c) ^ (b & c);
		std::uint32_t const temp2 = S0 + maj;
		hh = g; g = f; f = e; e = d + temp1;
		d = c; c = b; b = a; a = temp1 + temp2;
	}
	h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}

} // unnamed namespace

std::string Digest256::Hex() const
{
	static constexpr char digits[] = "0123456789abcdef";
	std::string out;
	out.reserve(64);
	for (auto b : bytes)
	{
		out.push_back(digits[b >> 4]);
		out.push_back(digits[b & 0xf]);
	}
	return out;
}

Digest256 Sha256(std::string_view data)
{
	// FIPS 180-4 initial hash: first 32 bits of the fractional parts of the square roots of
	// the first 8 primes.
	std::uint32_t h[8] = {
		0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
		0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
	};

	auto const* p = reinterpret_cast<std::uint8_t const*>(data.data());
	std::size_t remaining = data.size();
	while (remaining >= 64)
	{
		Compress(h, p);
		p += 64;
		remaining -= 64;
	}

	// Final one or two padded blocks: message tail, 0x80, zeros, 64-bit big-endian bit length.
	std::uint8_t block[64];
	std::memset(block, 0, sizeof(block));
	if (remaining)
		std::memcpy(block, p, remaining);
	block[remaining] = 0x80;
	if (remaining >= 56)
	{
		Compress(h, block);
		std::memset(block, 0, sizeof(block));
	}
	std::uint64_t const bit_length = std::uint64_t(data.size()) * 8;
	for (unsigned i = 0; i < 8; ++i)
		block[56 + i] = std::uint8_t(bit_length >> (56 - 8*i));
	Compress(h, block);

	Digest256 result;
	for (unsigned i = 0; i < 8; ++i)
	{
		result.bytes[4*i]   = std::uint8_t(h[i] >> 24);
		result.bytes[4*i+1] = std::uint8_t(h[i] >> 16);
		result.bytes[4*i+2] = std::uint8_t(h[i] >> 8);
		result.bytes[4*i+3] = std::uint8_t(h[i]);
	}
	return result;
}

} // namespace detail
} // namespace bertini
