//This file is part of Bertini 2.
//
//sha256_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//sha256_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with sha256_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Pins our clean-room SHA-256 to the FIPS 180-4 / NIST CAVP test vectors.  If any of
these ever fails, persistent System digests (ADR-0042) are broken -- do not "fix" the
expected strings; fix the implementation.
*/

#include <string>
#include "bertini2/detail/sha256.hpp"
#include <boost/test/unit_test.hpp>

using bertini::detail::Sha256;
using bertini::detail::Digest256;

BOOST_AUTO_TEST_SUITE(sha256)

// NIST FIPS 180-4 example: the empty message.
BOOST_AUTO_TEST_CASE(fips_empty_message)
{
	BOOST_CHECK_EQUAL(Sha256("").Hex(),
		"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

// NIST FIPS 180-4 example B.1: one-block message "abc".
BOOST_AUTO_TEST_CASE(fips_abc)
{
	BOOST_CHECK_EQUAL(Sha256("abc").Hex(),
		"ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}

// NIST FIPS 180-4 example B.2: two-block 448-bit message.
BOOST_AUTO_TEST_CASE(fips_two_block_message)
{
	BOOST_CHECK_EQUAL(
		Sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq").Hex(),
		"248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

// NIST CAVP long-message case: one million 'a's (exercises many blocks).
BOOST_AUTO_TEST_CASE(million_a)
{
	std::string const msg(1000000, 'a');
	BOOST_CHECK_EQUAL(Sha256(msg).Hex(),
		"cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

// Padding boundary cases: lengths 55, 56, 63, 64 straddle the one-vs-two final padded
// blocks decision (a 55-byte tail fits length in one block; 56+ forces a second).
BOOST_AUTO_TEST_CASE(padding_boundaries)
{
	// Expected values computed with coreutils sha256sum (independent implementation).
	BOOST_CHECK_EQUAL(Sha256(std::string(55, 'x')).Hex(),
		"d5e285683cd4efc02d021a5c62014694958901005d6f71e89e0989fac77e4072");
	BOOST_CHECK_EQUAL(Sha256(std::string(56, 'x')).Hex(),
		"04c26261370ee7541549d16dee320c723e3fd14671e66a099afe0a377c16888e");
	BOOST_CHECK_EQUAL(Sha256(std::string(63, 'x')).Hex(),
		"75220b47218278e656f2013bb8f0c455a25eaf01e86c64924e9d48d89776d6f2");
	BOOST_CHECK_EQUAL(Sha256(std::string(64, 'x')).Hex(),
		"7ce100971f64e7001e8fe5a51973ecdfe1ced42befe7ee8d5fd6219506b5393c");
}

BOOST_AUTO_TEST_CASE(digest_value_semantics)
{
	Digest256 const a = Sha256("abc");
	Digest256 const b = Sha256("abc");
	Digest256 const c = Sha256("abd");
	BOOST_CHECK(a == b);
	BOOST_CHECK(a != c);
	BOOST_CHECK((a < c) != (c < a));  // strict total order distinguishes them
	BOOST_CHECK_EQUAL(a.Hex().size(), 64u);
}

BOOST_AUTO_TEST_SUITE_END()
