//This file is part of Bertini 2.
//
//encoding_registry_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//encoding_registry_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with encoding_registry_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file The rule the encoding-version registries live by: AT MOST ONE ENCODING VERSION PER
RELEASE.  The two real registries are supposed to satisfy it, so checking them proves only
that a correct file passes.  These feed the rule registries that BREAK it, which is the half
that says the rule has teeth.

The rule exists because a version number that climbs with development rather than with
releases stops meaning anything: b2cfgenc/1 and /2 were both minted during 3.0.0's
development and neither ever shipped, so the number said "three encodings exist" when one
did.  Under the rule, an unreleased line is edited in place and the number counts encodings
that were really written.
*/

#include <string>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "test/utility/encoding_registry.hpp"

using bertini::testing::EncodingRegistryLine;
using bertini::testing::EncodingRegistryProblems;
using bertini::testing::ReleaseSeries;
using bertini::testing::kNeverShipped;

namespace {

// A registry line with one hash column, spelled the way the real files spell it.
EncodingRegistryLine Line(std::string version, std::string hash, std::string shipped)
{
    return EncodingRegistryLine{std::move(version), {std::move(hash)}, std::move(shipped)};
}

// The problems the rule finds, for a registry whose current version is its last line.
std::vector<std::string> Problems(std::vector<EncodingRegistryLine> const& lines,
                                  std::tuple<unsigned, unsigned, unsigned> release = {3, 5, 0})
{
    return EncodingRegistryProblems(lines, "b2enc/",
                                    lines.empty() ? std::string{} : lines.back().version, release);
}

} // unnamed namespace


BOOST_AUTO_TEST_SUITE(encoding_registry)

BOOST_AUTO_TEST_CASE(a_well_formed_registry_has_nothing_wrong_with_it)
{
    BOOST_CHECK(Problems({Line("b2enc/1", "aa", "3.0.0"),
                          Line("b2enc/2", "bb", "3.5.0")}).empty());
}


BOOST_AUTO_TEST_CASE(two_versions_in_one_release_is_the_thing_the_rule_forbids)
{
    // the failure mode the rule is named for: a development cycle minting version after version
    auto const problems = Problems({Line("b2enc/1", "aa", "3.0.0"),
                                    Line("b2enc/2", "bb", "3.5.0"),
                                    Line("b2enc/3", "cc", "3.5.0")});
    BOOST_REQUIRE_EQUAL(problems.size(), 1u);
    BOOST_CHECK(problems.front().find("same release") != std::string::npos);
}


BOOST_AUTO_TEST_CASE(releases_may_not_go_backwards)
{
    auto const problems = Problems({Line("b2enc/1", "aa", "3.5.0"),
                                    Line("b2enc/2", "bb", "3.0.0")});
    BOOST_CHECK_EQUAL(problems.size(), 1u);
}


BOOST_AUTO_TEST_CASE(a_version_may_not_claim_a_release_the_tree_has_not_reached)
{
    auto const problems = Problems({Line("b2enc/1", "aa", "3.6.0")}, {3, 5, 0});
    BOOST_REQUIRE_EQUAL(problems.size(), 1u);
    BOOST_CHECK(problems.front().find("has not reached") != std::string::npos);
}


BOOST_AUTO_TEST_CASE(a_never_shipped_version_may_only_precede_the_ones_that_shipped)
{
    // never-shipped versions are relics: before the rule, a dev cycle could supersede a version
    // before it reached a release.  Under the rule the unreleased line is edited instead, so a
    // new one cannot arise -- and one appearing AFTER a shipped version means the rule was broken.
    BOOST_CHECK(Problems({Line("b2enc/1", "aa", kNeverShipped),
                          Line("b2enc/2", "bb", "3.0.0"),
                          Line("b2enc/3", "cc", "3.5.0")}).empty());

    auto const problems = Problems({Line("b2enc/1", "aa", "3.0.0"),
                                    Line("b2enc/2", "bb", kNeverShipped)});
    BOOST_CHECK(!problems.empty());
}


BOOST_AUTO_TEST_CASE(version_numbers_are_dense_and_in_order)
{
    auto const skipped = Problems({Line("b2enc/1", "aa", "3.0.0"),
                                   Line("b2enc/3", "bb", "3.5.0")});
    BOOST_CHECK(!skipped.empty());

    auto const swapped = Problems({Line("b2enc/2", "aa", "3.0.0"),
                                   Line("b2enc/1", "bb", "3.5.0")});
    BOOST_CHECK(!swapped.empty());
}


BOOST_AUTO_TEST_CASE(the_last_line_must_be_what_this_build_writes)
{
    std::vector<EncodingRegistryLine> const lines{Line("b2enc/1", "aa", "3.0.0"),
                                                  Line("b2enc/2", "bb", "3.5.0")};
    BOOST_CHECK(EncodingRegistryProblems(lines, "b2enc/", "b2enc/2", {3, 5, 0}).empty());

    auto const stale = EncodingRegistryProblems(lines, "b2enc/", "b2enc/3", {3, 5, 0});
    BOOST_REQUIRE_EQUAL(stale.size(), 1u);
    BOOST_CHECK(stale.front().find("LAST registry line") != std::string::npos);
}


BOOST_AUTO_TEST_CASE(an_empty_registry_is_a_problem_not_a_pass)
{
    BOOST_CHECK_EQUAL(Problems({}).size(), 1u);
}


BOOST_AUTO_TEST_CASE(a_dev_version_is_the_release_it_is_heading_for)
{
    // 3.5.0.dev1 is the churn channel of 3.5.0, not a release of its own -- otherwise every dev
    // bump would look like a new release and the rule would permit a version per dev wheel
    BOOST_CHECK(ReleaseSeries("3.5.0.dev1") == ReleaseSeries("3.5.0"));
    BOOST_CHECK(ReleaseSeries("3.5.0rc1") == ReleaseSeries("3.5.0"));
    BOOST_CHECK(ReleaseSeries("3.4.0") < ReleaseSeries("3.5.0"));
    BOOST_CHECK(ReleaseSeries("3.5.0") < ReleaseSeries("3.5.1"));
}

BOOST_AUTO_TEST_SUITE_END()
