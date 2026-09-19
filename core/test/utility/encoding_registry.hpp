//This file is part of Bertini 2.
//
//test/utility/encoding_registry.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/utility/encoding_registry.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/utility/encoding_registry.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file test/utility/encoding_registry.hpp

\brief Reading and checking an encoding-version registry, shared by the config and system
identity tests (ADR-0043, ADR-0042).

A registry is the record of every canonical encoding the library has ever written, one line
per version, in order.  The rule the two tests share, and this header enforces:

**At most one encoding version per released library version.**  The version number is then the
count of encodings that have ever existed -- `b2sysenc/2` means "the second encoding b2 has ever
written" -- and a reader holding a record can tell which release wrote it.  Without the rule the
number drifts upward with development rather than with releases: `b2cfgenc/1` and `/2` were both
minted during 3.0.0's development and neither ever shipped, which is exactly the confusion the
rule prevents.

The practical consequence, and the reason this is not simply "append-only": while a release is
UNRELEASED, a further identity-affecting change updates the top line **in place** -- new hash,
same version token -- and regenerates the golden fixture.  Append-only still holds where it
matters, which is lines that shipped: those are history and are never edited.  The cost of the
in-place edit is that records written by an earlier build of the same unreleased line are not
recalled, which is the dev line behaving like a dev line.

Line format is `<version> <hash>... <first-shipped-in>`: the version token, one or more hash
columns (each file documents its own), and the release the version first shipped in.  Two
sentinels are allowed: `unrecorded` for a hash nobody captured at the time, and `never-shipped`
for a version that was superseded before any release used it.
*/

#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

namespace bertini {
namespace testing {

/// \brief The release sentinel for a version that was superseded before any release used it.
inline constexpr char kNeverShipped[] = "never-shipped";

/// \brief The hash sentinel for a value nobody captured when the line was written.
inline constexpr char kUnrecorded[] = "unrecorded";


/// \brief One line of an encoding-version registry.
struct EncodingRegistryLine
{
    std::string version;               ///< The version token, e.g. `b2cfgenc/4`.
    std::vector<std::string> hashes;   ///< The hash columns, in file order.
    std::string first_shipped_in;      ///< The release this version first shipped in, or `never-shipped`.
};


/**
\brief Read a registry file.

\param path The registry file.
\param num_hash_columns How many hash columns this registry carries between the version and the
       release (the config registry has two -- encoders and composition; the system registry has one).
\return The lines, in file order.
*/
inline std::vector<EncodingRegistryLine> ReadEncodingRegistry(std::filesystem::path const& path,
                                                              std::size_t num_hash_columns)
{
    std::ifstream in(path);
    if (!in.good())
        throw std::runtime_error("version registry missing: " + path.string());

    std::vector<EncodingRegistryLine> lines;
    std::string raw;
    for (unsigned line_number = 1; std::getline(in, raw); ++line_number)
    {
        if (raw.empty() || raw[0] == '#')
            continue;
        std::istringstream fields(raw);
        EncodingRegistryLine line;
        if (!(fields >> line.version))
            continue;
        line.hashes.resize(num_hash_columns);
        for (auto& h : line.hashes)
            if (!(fields >> h))
                throw std::runtime_error(path.string() + ":" + std::to_string(line_number)
                    + ": expected " + std::to_string(num_hash_columns)
                    + " hash column(s) then the release it first shipped in");
        if (!(fields >> line.first_shipped_in))
            throw std::runtime_error(path.string() + ":" + std::to_string(line_number)
                + ": no release column -- every line says which release it first shipped in, or '"
                + kNeverShipped + "'");
        lines.push_back(std::move(line));
    }
    return lines;
}


/// \brief A release as a comparable (major, minor, patch); pre-release suffixes are dropped, so
/// `3.5.0.dev1` and `3.5.0` are the same release series -- a dev wheel is the churn channel of
/// the release it is heading for, not a release of its own.
inline std::tuple<unsigned, unsigned, unsigned> ReleaseSeries(std::string const& text)
{
    unsigned part[3] = {0, 0, 0};
    std::size_t at = 0;
    for (unsigned ii = 0; ii < 3 && at < text.size(); ++ii)
    {
        std::size_t digits = 0;
        while (at + digits < text.size() && std::isdigit(static_cast<unsigned char>(text[at + digits])))
            ++digits;
        if (digits == 0)
            break;
        part[ii] = static_cast<unsigned>(std::stoul(text.substr(at, digits)));
        at += digits;
        if (at < text.size() && text[at] == '.')
            ++at;
        else
            break;
    }
    return {part[0], part[1], part[2]};
}


/// \brief The release series this working tree is building, from the top-level VERSION file.
inline std::tuple<unsigned, unsigned, unsigned> CurrentReleaseSeries()
{
    auto const version_path =
        std::filesystem::path(__FILE__).parent_path() / ".." / ".." / ".." / "VERSION";
    std::ifstream in(version_path);
    if (!in.good())
        throw std::runtime_error("cannot read the VERSION file at " + version_path.string());
    std::string text;
    std::getline(in, text);
    return ReleaseSeries(text);
}


/**
\brief Everything wrong with a registry that does not depend on which encoding it describes.

Dense sequential version numbers; the last line is the version this build writes; at most one
line per release, in release order; and nothing claiming a release the tree has not reached.
The hash columns are the caller's business -- each test knows how to recompute its own.

A pure function rather than a wall of assertions, so the rule itself can be tested against
registries that break it (encoding_registry_test.cpp) rather than only against the two real
files, which are supposed to be correct.

\param lines The registry, as read.
\param version_prefix The token prefix, e.g. `b2cfgenc/`.
\param current_version The version token this build writes.
\param current_release The release series the tree is building.
\return One sentence per problem; empty when the registry is well formed.
*/
inline std::vector<std::string> EncodingRegistryProblems(
    std::vector<EncodingRegistryLine> const& lines,
    std::string const& version_prefix,
    std::string const& current_version,
    std::tuple<unsigned, unsigned, unsigned> const& current_release)
{
    std::vector<std::string> problems;
    if (lines.empty())
    {
        problems.push_back("the registry is empty");
        return problems;
    }

    for (std::size_t ii = 0; ii < lines.size(); ++ii)
        if (lines[ii].version != version_prefix + std::to_string(ii + 1))
            problems.push_back("line " + std::to_string(ii + 1) + " is " + lines[ii].version
                + ", but versions are dense and in order, so it must be "
                + version_prefix + std::to_string(ii + 1));

    if (lines.back().version != current_version)
        problems.push_back("the LAST registry line must be the version this build writes ("
            + current_version + "); found " + lines.back().version);

    // one line per release, in order.  A never-shipped version is a relic of the days before
    // this rule; a new one cannot arise, because an unreleased line is edited rather than
    // superseded, so they may only appear before the first line that did ship.
    bool seen_shipped = false;
    std::tuple<unsigned, unsigned, unsigned> previous{0, 0, 0};
    std::string previous_text;
    for (auto const& line : lines)
    {
        if (line.first_shipped_in == kNeverShipped)
        {
            if (seen_shipped)
                problems.push_back(line.version + " is marked " + kNeverShipped
                    + " but an earlier version already shipped -- under one-version-per-release "
                      "an unreleased version is edited in place, never superseded, so a "
                      "never-shipped version cannot follow a released one");
            continue;
        }
        seen_shipped = true;
        auto const series = ReleaseSeries(line.first_shipped_in);
        if (!previous_text.empty() && !(series > previous))
            problems.push_back("two encoding versions claim the same release, or go backwards: "
                + line.version + " says " + line.first_shipped_in + " after " + previous_text
                + ".  At most ONE encoding version per release -- while a release is unreleased, "
                  "an identity-affecting change updates the LAST line in place (new hash, same "
                  "token) and regenerates the golden fixture, rather than minting another version");
        previous = series;
        previous_text = line.first_shipped_in;
    }

    if (lines.back().first_shipped_in != kNeverShipped
        && ReleaseSeries(lines.back().first_shipped_in) > current_release)
        problems.push_back(lines.back().version + " claims to have first shipped in "
            + lines.back().first_shipped_in + ", which this working tree has not reached");

    return problems;
}


/// \brief Report every problem EncodingRegistryProblems finds with one of the real registries.
/// \param lines The registry, as read.
/// \param version_prefix The token prefix, e.g. `b2cfgenc/`.
/// \param current_version The version token this build writes.
/// \param path The registry path, for messages.
inline void CheckEncodingRegistryShape(std::vector<EncodingRegistryLine> const& lines,
                                       std::string const& version_prefix,
                                       std::string const& current_version,
                                       std::filesystem::path const& path)
{
    auto const problems = EncodingRegistryProblems(lines, version_prefix, current_version,
                                                   CurrentReleaseSeries());
    for (auto const& p : problems)
        BOOST_ERROR(path.string() + ": " + p);
    BOOST_REQUIRE_MESSAGE(problems.empty() && !lines.empty(),
        "the registry is not well formed, so its hash columns cannot be trusted: " << path);
}


/**
\brief The advice a failing hash column should give, which depends on whether the top line shipped.

An unreleased top line is edited in place; a shipped one must be superseded by a new version.

\param lines The registry, as read.
\param version_prefix The token prefix, e.g. `b2cfgenc/`.
\param column A human name for the hash column that moved.
\param found The hash this build computes.
\return A sentence naming the file edit to make.
*/
inline std::string AdviceForMovedHash(std::vector<EncodingRegistryLine> const& lines,
                                      std::string const& version_prefix,
                                      std::string const& column,
                                      std::string const& found)
{
    bool const top_line_shipped = !lines.empty()
        && lines.back().first_shipped_in != kNeverShipped
        && ReleaseSeries(lines.back().first_shipped_in) <= CurrentReleaseSeries()
        && ReleaseSeries(lines.back().first_shipped_in) != CurrentReleaseSeries();

    if (top_line_shipped)
        return "the last registry line belongs to a RELEASED version, so it is history and is not "
               "edited: bump the version token and APPEND a line {"
               + version_prefix + std::to_string(lines.size() + 1) + ", " + found
               + ", <this release>}, and regenerate the golden fixture, all in the same commit.";

    return "the last registry line belongs to this UNRELEASED release, so it is the line to "
           "change: replace its " + column + " with " + found
           + " -- same version token -- and regenerate the golden fixture, in the same commit.";
}


} // namespace testing
} // namespace bertini
