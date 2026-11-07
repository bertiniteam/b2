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
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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
System HomotopyFromTheCurrentState();
System SeededHomotopy(unsigned long seed)
{
    SetGlobalSeed(seed);
    return HomotopyFromTheCurrentState();
}

// ... drawing from wherever the random state is now
System HomotopyFromTheCurrentState()
{
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

// The effective per-solve seed: AdvanceSolveSeed takes the next link of a chain rooted at
// the master seed, so a seedless solve records a seed that reproduces it STANDALONE instead
// of a session master that silently under-determines a mid-session run.
BOOST_AUTO_TEST_CASE(solve_seeds_advance_along_a_chain_from_the_master_seed)
{
    SetGlobalSeed(42);
    auto const s1 = AdvanceSolveSeed();
    auto const s2 = AdvanceSolveSeed();
    BOOST_CHECK(s1 != 0u);                       // 0 means "entropy" to a classic input file
    BOOST_CHECK(s1 != s2);                       // consecutive solves get distinct seeds
    BOOST_CHECK(s1 <= 0xFFFFFFFFul);             // 32-bit portable (Windows unsigned long)
    BOOST_CHECK(s2 <= 0xFFFFFFFFul);

    SetGlobalSeed(42);                           // same master =>
    BOOST_CHECK_EQUAL(AdvanceSolveSeed(), s1);   // the same seeds, in order
    BOOST_CHECK_EQUAL(AdvanceSolveSeed(), s2);

    // DELIBERATELY independent of the draw streams: a multithreaded solve leaves the
    // calling thread's stream at a scheduling-dependent position, so the derivation
    // must not read it -- intervening draws change nothing
    SetGlobalSeed(42);
    (void)RandomRat();
    BOOST_CHECK_EQUAL(AdvanceSolveSeed(), s1);

    SetGlobalSeed(43);                           // different master, different chain
    BOOST_CHECK(AdvanceSolveSeed() != s1);
}

// the former name advances the same chain until it is removed in 5.0
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
BOOST_AUTO_TEST_CASE(the_former_name_derive_solve_seed_advances_the_same_chain)
{
    SetGlobalSeed(42);
    auto const first = AdvanceSolveSeed();
    auto const second = AdvanceSolveSeed();
    SetGlobalSeed(42);
    BOOST_CHECK_EQUAL(DeriveSolveSeed(), first);
    BOOST_CHECK_EQUAL(AdvanceSolveSeed(), second);
}
#pragma GCC diagnostic pop

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

// ---- seed 0 (#475, ADR-0072) ----

// 0 is a seed like any other, as in numpy and Python's random.  It used to mean "draw one from
// entropy", so a script that set seed 0 to make a run reproducible got a different run each time.
BOOST_AUTO_TEST_CASE(seed_zero_is_an_ordinary_seed)
{
    SetGlobalSeed(0);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), 0ul);
    auto const r = RandomRat();
    auto const m = RandomMp(30);
    auto const digest = SeededHomotopy(0).ContentDigest().Hex();

    SetGlobalSeed(0);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), 0ul);
    BOOST_CHECK_EQUAL(RandomRat(), r);
    BOOST_CHECK_EQUAL(RandomMp(30), m);
    BOOST_CHECK_EQUAL(SeededHomotopy(0).ContentDigest().Hex(), digest);

    BOOST_CHECK(SeededHomotopy(1).ContentDigest().Hex() != digest);
}

// Entropy is its own request, and it says which seed it chose, so the run can be reproduced.
BOOST_AUTO_TEST_CASE(a_seed_from_entropy_is_reported_and_reproduces_the_run)
{
    auto const chosen = SetGlobalSeedFromEntropy();
    BOOST_CHECK_EQUAL(GetGlobalSeed(), chosen);
    auto const r = RandomRat();

    SetGlobalSeed(chosen);
    BOOST_CHECK_EQUAL(RandomRat(), r);
}

// The classic input format keeps Bertini 1's meaning: randomseed 0, its default, draws a seed
// from entropy; any other value is the seed.
BOOST_AUTO_TEST_CASE(classic_randomseed_zero_draws_from_entropy)
{
    BOOST_CHECK_EQUAL(ApplyClassicRandomSeed(5), 5ul);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), 5ul);

    auto const first = ApplyClassicRandomSeed(0);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), first);
    auto const second = ApplyClassicRandomSeed(0);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), second);
    // two draws from entropy agree with probability about 2^-32
    BOOST_CHECK(first != second);
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

    // re-pinned for b2sysenc/3 (the content-key tie-break, ADR-0073), which is edited in place
    // while 4.1.0 is unreleased (ADR-0061), so this moves whenever the encoding does until release
    char const* const pinned = "2508128366cab73dd135a084eeba2a28dec4866bf11205874d32c6fc6e5a1a0a";
    if (hex != pinned)
    {
        std::cout << "[diagnostic] canonical encoding of the mismatching homotopy follows\n"
                  << "-----8<-----\n" << homotopy.CanonicalEncodingText() << "\n-----8<-----\n";
    }
    BOOST_CHECK_EQUAL(hex, pinned);
}

// ---- getting and setting the random state (bertiniteam/b2#518) ----

namespace {

// a mix of every kind of draw: multiprecision numbers, a double, the stream's raw bytes
// (21 bits, so the position ends partway into a block), and the solve-seed chain
std::string DrawsFromHere()
{
    std::ostringstream out;
    out << RandomMp(30) << ' ' << records::ThreadDrawStream().Bits(21) << ' '
        << RandReal() << ' ' << AdvanceSolveSeed() << ' ' << RandomMp(60);
    return out.str();
}

} // unnamed namespace

BOOST_AUTO_TEST_CASE(a_state_set_again_draws_what_followed_getting_it)
{
    SetGlobalSeed(1042);
    DrawsFromHere();
    auto const state = GetRandomState();
    auto const after_getting_it = DrawsFromHere();

    SetGlobalSeed(77);              // another seed, another solve ordinal, other positions
    DrawsFromHere();
    SetRandomState(state);
    BOOST_CHECK_EQUAL(DrawsFromHere(), after_getting_it);
    BOOST_CHECK_EQUAL(GetGlobalSeed(), 1042ul);
}

BOOST_AUTO_TEST_CASE(a_homotopy_built_after_setting_a_state_is_the_one_built_after_getting_it)
{
    SetGlobalSeed(42);
    DrawsFromHere();
    auto const state = GetRandomState();
    auto const after_getting_it = HomotopyFromTheCurrentState().ContentDigest().Hex();

    HomotopyFromTheCurrentState();
    SetRandomState(state);
    BOOST_CHECK_EQUAL(HomotopyFromTheCurrentState().ContentDigest().Hex(), after_getting_it);
}

BOOST_AUTO_TEST_CASE(a_state_holds_the_seed_and_ordinal_and_compares_by_its_fields)
{
    SetGlobalSeed(5);
    DrawsFromHere();
    auto const state = GetRandomState();
    BOOST_CHECK_EQUAL(state.global_seed, 5ul);
    BOOST_CHECK(state.seed_chosen);
    BOOST_CHECK_EQUAL(state.solve_ordinal, 1u);
    BOOST_CHECK(GetRandomState() == state);

    DrawsFromHere();
    BOOST_CHECK(GetRandomState() != state);
    auto moved_stream = state;
    moved_stream.draw_stream.Uint64();
    BOOST_CHECK(moved_stream != state);
}

BOOST_AUTO_TEST_CASE(a_state_built_field_by_field_draws_as_one_got)
{
    SetGlobalSeed(11);
    DrawsFromHere();
    auto const got = GetRandomState();
    auto const from_got = DrawsFromHere();

    RandomState built;
    built.global_seed = got.global_seed;
    built.seed_chosen = got.seed_chosen;
    built.solve_ordinal = got.solve_ordinal;
    built.draw_stream = got.draw_stream;
    SetGlobalSeed(12);
    SetRandomState(built);
    BOOST_CHECK_EQUAL(DrawsFromHere(), from_got);
}

BOOST_AUTO_TEST_CASE(a_state_written_as_text_reads_back_equal)
{
    SetGlobalSeed(5);
    DrawsFromHere();
    auto const state = GetRandomState();
    auto const text = state.ToText();
    BOOST_CHECK_EQUAL(text.rfind("b2randstate/1\nglobal_seed 5\nseed_chosen 1\nsolve_ordinal 1\n"
                                 "draw_stream 1 ", 0), 0u);
    BOOST_CHECK(RandomState::FromText(text) == state);
    BOOST_CHECK_EQUAL(RandomState::FromText(text).ToText(), text);
}

BOOST_AUTO_TEST_CASE(a_text_that_is_not_a_state_is_refused)
{
    SetGlobalSeed(9);
    DrawsFromHere();
    auto const text = GetRandomState().ToText();

    std::string other_version = text;
    other_version.replace(0, std::string("b2randstate/1").size(), "b2randstate/2");
    std::string short_key = text;
    auto const at = short_key.find("draw_stream 1 ") + std::string("draw_stream 1 ").size();
    short_key.erase(at, 2);
    std::string no_stream = text.substr(0, text.find("draw_stream"));

    for (auto const& bad : {std::string(), std::string("global_seed 9"), other_version, short_key,
                            no_stream, text + "more\n"})
        BOOST_CHECK_THROW(RandomState::FromText(bad), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(a_state_goes_through_text_and_binary_archives)
{
    SetGlobalSeed(13);
    DrawsFromHere();
    auto const state = GetRandomState();

    std::stringstream text_archive;
    {
        boost::archive::text_oarchive out(text_archive);
        out << state;
    }
    RandomState from_text;
    {
        boost::archive::text_iarchive in(text_archive);
        in >> from_text;
    }
    BOOST_CHECK(from_text == state);

    std::stringstream binary_archive;
    {
        boost::archive::binary_oarchive out(binary_archive);
        out << state;
    }
    RandomState from_binary;
    {
        boost::archive::binary_iarchive in(binary_archive);
        in >> from_binary;
    }
    BOOST_CHECK(from_binary == state);
}

BOOST_AUTO_TEST_SUITE_END()
