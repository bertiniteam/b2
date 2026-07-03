//This file is part of Bertini 2.
//
//config_encoding.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//config_encoding.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with config_encoding.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file config_encoding.hpp

\brief Canonical exact encoding + persistent digests for configuration structs
(ADR-0043; the structured-output-directory arc, rung 1).

A solver's settings are half of a computation's ask identity: records reference the
configuration a path was tracked under by digest, so equal settings must digest equally
FOREVER -- across runs, compilers, machines, and versions.  Hence the same rules as the
System encoding (ADR-0042):

- versioned header (`b2cfgenc/1`) baked into every digest -- a spec change is a new
  keyspace, never silent drift;
- exact values only: doubles encode as their IEEE-754 bit pattern (`d64:<16 hex>`),
  NEVER decimal round-trips; rationals via exact `.str()`; enums via fixed string
  tables (never numeric values or typeid); strings as netstrings;
- one encoder per struct, listing fields by name in declaration order.  Adding a field
  to a config struct REQUIRES extending its encoder (and bumping `b2cfgenc` if the
  change is identity-affecting) -- the golden fixture in config_digest_test.cpp fails
  loudly otherwise.

Deliberately excluded from encoding: `algorithm::RandomConfig` (the seed is its own slot
in the ask identity, beside the config digest -- never inside it) and
`ZeroDimConfig::num_threads` (thread count must not change what was computed).
*/

#pragma once

#include <string>

#include "bertini2/detail/sha256.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"

namespace bertini {
namespace records {

/// \brief The version tag baked into every config digest; bump on any encoding change.
constexpr char ConfigEncodingVersion[] = "b2cfgenc/1";

/// \brief The fixed canonical name of a predictor choice (never its numeric value).
std::string CanonicalName(tracking::Predictor p);
/// \brief The fixed canonical name of an algorithm (track-type) choice.
std::string CanonicalName(algorithm::classic::AlgoChoice a);
/// \brief The fixed canonical name of an endgame choice.
std::string CanonicalName(algorithm::classic::EndgameChoice e);

/// \brief Canonical encoding of stepping settings.
std::string CanonicalEncoding(tracking::SteppingConfig const& c);
/// \brief Canonical encoding of Newton-corrector settings.
std::string CanonicalEncoding(tracking::NewtonConfig const& c);
/// \brief Canonical encoding of fixed-precision settings.
std::string CanonicalEncoding(tracking::FixedPrecisionConfig const& c);
/// \brief Canonical encoding of adaptive-multiple-precision settings.
std::string CanonicalEncoding(tracking::AdaptiveMultiplePrecisionConfig const& c);
/// \brief Canonical encoding of a predictor choice (as a one-field pseudo-config).
std::string CanonicalEncoding(tracking::Predictor p);

/// \brief Canonical encoding of endgame security (divergence-bailout) settings.
std::string CanonicalEncoding(endgame::SecurityConfig const& c);
/// \brief Canonical encoding of common endgame settings.
std::string CanonicalEncoding(endgame::EndgameConfig const& c);
/// \brief Canonical encoding of power-series-endgame settings.
std::string CanonicalEncoding(endgame::PowerSeriesConfig const& c);
/// \brief Canonical encoding of Cauchy-endgame settings.
std::string CanonicalEncoding(endgame::CauchyConfig const& c);
/// \brief Canonical encoding of track-back endgame settings.
std::string CanonicalEncoding(endgame::TrackBackConfig const& c);

/// \brief Canonical encoding of Newton/path tolerances.
std::string CanonicalEncoding(algorithm::TolerancesConfig const& c);
/// \brief Canonical encoding of midpath (path-crossing check) settings.
std::string CanonicalEncoding(algorithm::MidPathConfig const& c);
/// \brief Canonical encoding of crossed-path auto-retrack settings.
std::string CanonicalEncoding(algorithm::AutoRetrackConfig const& c);
/// \brief Canonical encoding of solution-sharpening settings.
std::string CanonicalEncoding(algorithm::SharpeningConfig const& c);
/// \brief Canonical encoding of regeneration settings.
std::string CanonicalEncoding(algorithm::RegenerationConfig const& c);
/// \brief Canonical encoding of post-processing (classification) settings.
std::string CanonicalEncoding(algorithm::PostProcessingConfig const& c);
/// \brief Canonical encoding of top-level zero-dim solve settings.  Excludes
/// num_threads (transient: thread count must not change identity).
std::string CanonicalEncoding(algorithm::ZeroDimConfig const& c);
/// \brief Canonical encoding of the algorithm (track-type) selection.
std::string CanonicalEncoding(algorithm::MetaConfig const& c);
/// \brief Canonical encoding of the endgame selection.
std::string CanonicalEncoding(algorithm::classic::EndgameChoiceConfig const& c);

/**
\brief The persistent digest of one configuration: SHA-256 over its versioned canonical
encoding.

\tparam ConfigT Any type with a CanonicalEncoding overload above.
\param config The configuration to digest.
\return The 256-bit digest (stable cross-session; the config half of an ask identity).
*/
template <typename ConfigT>
detail::Digest256 ConfigDigest(ConfigT const& config)
{
	return detail::Sha256(std::string(ConfigEncodingVersion) + "\n" + CanonicalEncoding(config));
}

/**
\brief One digest for a solver's FULL settings: the versioned concatenation of several
configs' canonical encodings, in the order given.

Order matters and is part of the contract: a given solver kind should always compose
its settings digest in one fixed order (document it at the call site).

\tparam ConfigTs Types each having a CanonicalEncoding overload.
\param configs The configurations, in the solver's fixed order.
\return The 256-bit digest of the combined settings.
*/
template <typename... ConfigTs>
detail::Digest256 SettingsDigest(ConfigTs const&... configs)
{
	std::string combined(ConfigEncodingVersion);
	((combined += "\n" + CanonicalEncoding(configs)), ...);
	return detail::Sha256(combined);
}

/**
\brief Render one or more canonical config encodings as pretty JSON, for archiving in a
structured output directory.

The DIGEST contract stays on the canonical text (`b2cfgenc/<n>`); this JSON view is
derived from it mechanically -- exact values preserved as strings -- because JSON is the
format humans and tools in the directory already speak.

\param canonical_text Newline-separated canonical encodings, first line the version tag.
\param digest_hex The settings digest these encodings produce (recorded alongside).
\return A pretty-printed JSON document: {"schema", "digest", "configs": {Name: {k: v}}}.
*/
std::string ConfigTextAsJson(std::string const& canonical_text, std::string const& digest_hex);

} // namespace records
} // namespace bertini
