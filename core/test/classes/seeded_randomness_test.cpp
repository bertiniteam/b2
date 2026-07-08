//This file is part of Bertini 2.
//
//seeded_randomness_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//seeded_randomness_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with seeded_randomness_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for seed-rooted randomness (b2rand/1, ADR-0044).  The literal expected
values here ARE the derivation contract: if one fails after a draw-layer change, that
change redefines what every seed means -- bump b2rand/<n>, never edit the expectation.
Because they are checked on all three CI platforms, they also enforce cross-platform
draw reproducibility.
*/

#include <iostream>
#include <sstream>

#include <boost/test/unit_test.hpp>

#include "bertini2/records/derive.hpp"
#include "bertini2/random.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

using namespace bertini;

namespace {

System ParseCircleLine()
{
	std::string const str = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";
	System sys;
	[[maybe_unused]] bool success = parsing::classic::parse(str.begin(), str.end(), sys);
	return sys;
}

// build the whole random cascade for one seed: TD start + gamma homotopy + patch
System SeededHomotopy(unsigned long seed)
{
	SetGlobalSeed(seed);
	auto sys = ParseCircleLine();
	sys.Homogenize();
	sys.AutoPatch();
	start_system::TotalDegreeBinomial start(sys);
	return MakeHomotopy(sys, start, "t");   // null gamma -> a seeded random gamma
}

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(seeded_randomness)

// ---- the primitive contract (b2rand/1): pinned raw draws ----

BOOST_AUTO_TEST_CASE(stream_draws_are_the_pinned_contract)
{
	records::DrawStream stream;
	stream.Reseed(42, 0, 0);

	std::ostringstream first;
	first << std::hex << stream.Uint64();
	BOOST_CHECK_EQUAL(first.str(), "32ea87f405f7d5de");

	stream.Reseed(42, 0, 0);
	(void)stream.Uint64();
	std::ostringstream unit;
	unit.precision(17);
	unit << stream.UnitDouble();
	BOOST_CHECK_EQUAL(unit.str(), "0.64174190966992761");
}

BOOST_AUTO_TEST_CASE(seeded_session_draws_are_pinned)
{
	SetGlobalSeed(42);
	BOOST_CHECK_EQUAL(RandomRat().str(), "14324747069416753937984000168135157906230970132321/36331215644497613286232963764556423811994956822371");
	SetGlobalSeed(42);
	BOOST_CHECK_EQUAL(RandomMp(30).str(), "0.32657335689626739904012359602498");
}

// ---- reproducibility semantics ----

// The effective per-solve seed: DeriveSolveSeed captures the session stream's current
// position, so a seedless solve records a seed that reproduces it STANDALONE instead
// of a session master that silently under-determines a mid-session run.
BOOST_AUTO_TEST_CASE(derived_solve_seeds_capture_the_stream_deterministically)
{
	SetGlobalSeed(42);
	auto const s1 = DeriveSolveSeed();
	auto const s2 = DeriveSolveSeed();
	BOOST_CHECK(s1 != 0u);                       // 0 means "entropy" to SetGlobalSeed
	BOOST_CHECK(s1 != s2);                       // consecutive solves get distinct seeds
	BOOST_CHECK(s1 <= 0xFFFFFFFFul);             // 32-bit portable (Windows unsigned long)
	BOOST_CHECK(s2 <= 0xFFFFFFFFul);

	SetGlobalSeed(42);                           // same master =>
	BOOST_CHECK_EQUAL(DeriveSolveSeed(), s1);    // the same derived seeds, in order
	BOOST_CHECK_EQUAL(DeriveSolveSeed(), s2);

	// DELIBERATELY independent of the draw streams: a multithreaded solve leaves the
	// calling thread's stream at a scheduling-dependent position, so the derivation
	// must not read it -- intervening draws change nothing
	SetGlobalSeed(42);
	(void)RandomRat();
	BOOST_CHECK_EQUAL(DeriveSolveSeed(), s1);

	SetGlobalSeed(43);                           // different master, different chain
	BOOST_CHECK(DeriveSolveSeed() != s1);
}

BOOST_AUTO_TEST_CASE(same_seed_same_draws_different_seed_different_draws)
{
	SetGlobalSeed(42);
	auto const a1 = RandomRat();
	auto const a2 = RandomMp(30);

	SetGlobalSeed(42);
	BOOST_CHECK_EQUAL(RandomRat(), a1);
	BOOST_CHECK_EQUAL(RandomMp(30), a2);

	SetGlobalSeed(43);
	BOOST_CHECK(RandomRat() != a1);
}

BOOST_AUTO_TEST_CASE(per_path_streams_are_deterministic_and_domain_separated)
{
	SetGlobalSeed(42);
	ReseedThisThread(7);
	auto const path7_draw = RandomMp(30);

	ReseedThisThread(7);
	BOOST_CHECK_EQUAL(RandomMp(30), path7_draw);         // same path, same draws

	ReseedThisThread(8);
	BOOST_CHECK(RandomMp(30) != path7_draw);             // different path stream

	SetGlobalSeed(42);                                    // the setup stream
	BOOST_CHECK(RandomMp(30) != path7_draw);             // never collides with a path stream
}

// ---- issue #294: the friendly factory draws are seeded and real-when-asked ----

// The bounded-modulus factory draws (behind bertini.random_real / random_complex / random_vector)
// are continuous, seed-reproducible, and -- for the real draw -- genuinely real.  This is the
// generic-direction tool the notebook wanted, distinct from the quantized orthonormal random_matrix.
BOOST_AUTO_TEST_CASE(bounded_modulus_factories_are_seeded_and_real_when_asked)
{
	using bertini::multiprecision::RandomRealBoundedModulus;
	using bertini::multiprecision::RandomComplexBoundedModulus;

	SetGlobalSeed(42);
	auto const r1 = RandomRealBoundedModulus();
	auto const c1 = RandomComplexBoundedModulus();

	BOOST_CHECK_EQUAL(r1.imag(), real_mp(0));            // a real draw is actually real

	SetGlobalSeed(42);                                   // same seed reproduces the draws, in order
	BOOST_CHECK_EQUAL(RandomRealBoundedModulus(), r1);
	BOOST_CHECK_EQUAL(RandomComplexBoundedModulus(), c1);

	SetGlobalSeed(43);                                   // a different seed => a different draw
	BOOST_CHECK(RandomRealBoundedModulus() != r1);       // (not the seed-independent orthonormal footgun)
}

// ---- the acceptance test: same seed => digest-identical homotopies ----

BOOST_AUTO_TEST_CASE(same_seed_builds_digest_identical_homotopies)
{
	auto const first = SeededHomotopy(42).ContentDigest();
	auto const again = SeededHomotopy(42).ContentDigest();
	BOOST_CHECK_EQUAL(first.Hex(), again.Hex());

	auto const other = SeededHomotopy(43).ContentDigest();
	BOOST_CHECK(first.Hex() != other.Hex());
}

// The cross-platform oracle: this digest folds every seeded draw in the construction
// cascade (TD coefficients, gamma, patch) through mpfr arithmetic; CI checks it on
// Linux, macOS, and Windows, so a platform-dependent draw or rounding difference
// fails loudly here.
BOOST_AUTO_TEST_CASE(seed_42_homotopy_digest_is_pinned_cross_platform)
{
	// Diagnostic stage digests (CI, 2026-07-03): the pinned value was minted on arm64 and
	// x86_64 CI computes a DIFFERENT digest (6e7155e9...), while the raw-draw pins and the
	// parsed-system fixtures pass everywhere -- so some stage of the construction cascade
	// is architecture-dependent.  Print each stage's digest, and on final mismatch dump the
	// full canonical encoding, so the CI log pinpoints the diverging stage and bytes.
	SetGlobalSeed(42);
	auto sys = ParseCircleLine();
	std::cout << "[stage] parsed:      " << sys.ContentDigest().Hex() << "\n";
	sys.Homogenize();
	std::cout << "[stage] homogenized: " << sys.ContentDigest().Hex() << "\n";
	sys.AutoPatch();
	std::cout << "[stage] patched:     " << sys.ContentDigest().Hex() << "\n";
	start_system::TotalDegreeBinomial start(sys);
	std::cout << "[stage] td start:    " << start.ContentDigest().Hex() << "\n";
	auto homotopy = MakeHomotopy(sys, start, "t");
	auto const hex = homotopy.ContentDigest().Hex();
	std::cout << "[stage] homotopy:    " << hex << "\n";

	char const* const pinned = "4b2970766e8b58fb50f03fe9d7ce5bec58c07d5897b8637e55838248fc4fd2ac";
	if (hex != pinned)
	{
		std::cout << "[diagnostic] canonical encoding of the mismatching homotopy follows\n"
		          << "-----8<-----\n" << homotopy.CanonicalEncodingText() << "\n-----8<-----\n";
	}
	BOOST_CHECK_EQUAL(hex, pinned);
}

BOOST_AUTO_TEST_SUITE_END()
