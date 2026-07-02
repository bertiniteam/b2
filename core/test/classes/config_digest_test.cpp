//This file is part of Bertini 2.
//
//config_digest_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//config_digest_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with config_digest_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for config canonical encoding + digests (ADR-0043) and the CROSS-SESSION
GOLDEN FIXTURE (data/config_digest_fixture.txt).  If a golden digest changes, the
config encoding drifted: bump `b2cfgenc/<n>` and regenerate the fixture in the same
commit, never silently.
*/

#include <filesystem>
#include <fstream>
#include <map>
#include <string>

#include <boost/test/unit_test.hpp>

#include "bertini2/records/config_encoding.hpp"

using namespace bertini;
using records::CanonicalEncoding;
using records::ConfigDigest;
using records::SettingsDigest;

BOOST_AUTO_TEST_SUITE(config_digests)

// ---- semantics ----

BOOST_AUTO_TEST_CASE(equal_configs_digest_equally_across_independent_construction)
{
	tracking::SteppingConfig a, b;
	BOOST_CHECK(ConfigDigest(a) == ConfigDigest(b));

	endgame::CauchyConfig ca, cb;
	BOOST_CHECK(ConfigDigest(ca) == ConfigDigest(cb));
}

BOOST_AUTO_TEST_CASE(every_field_change_moves_the_digest)
{
	// stepping: touch each field in turn; the digest must move each time
	tracking::SteppingConfig base;
	auto const base_digest = ConfigDigest(base);

	auto touched = base; touched.initial_step_size = mpq_rational(1, 7);
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.max_step_size = mpq_rational(1, 7);
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.min_step_size = 1e-99;
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.step_size_success_factor = mpq_rational(3, 1);
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.step_size_fail_factor = mpq_rational(1, 3);
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.consecutive_successful_steps_before_stepsize_increase = 6;
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.min_num_steps = 2;
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.max_num_steps = 99;
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
	touched = base; touched.frequency_of_CN_estimation = 7;
	BOOST_CHECK(!(ConfigDigest(touched) == base_digest));
}

BOOST_AUTO_TEST_CASE(doubles_encode_exactly_not_by_decimal)
{
	// two doubles that print identically at low decimal precision must still differ
	algorithm::TolerancesConfig a, b;
	a.final_tolerance = 1e-11;
	b.final_tolerance = 1e-11 * (1 + 1e-16);   // one ulp-ish away
	if (a.final_tolerance != b.final_tolerance)  // guard: only if genuinely distinct doubles
		BOOST_CHECK(!(ConfigDigest(a) == ConfigDigest(b)));

	// and the encoding is the bit pattern, not a formatted number
	BOOST_CHECK(CanonicalEncoding(a).find("d64:") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(num_threads_is_not_identity)
{
	algorithm::ZeroDimConfig a, b;
	a.num_threads = 1;
	b.num_threads = 64;
	BOOST_CHECK(ConfigDigest(a) == ConfigDigest(b));   // transient: excluded by design
}

BOOST_AUTO_TEST_CASE(settings_digest_composes_and_order_matters)
{
	tracking::SteppingConfig s;
	tracking::NewtonConfig n;
	BOOST_CHECK(SettingsDigest(s, n) == SettingsDigest(s, n));
	BOOST_CHECK(!(SettingsDigest(s, n) == SettingsDigest(n, s)));  // order is contract
	BOOST_CHECK(!(SettingsDigest(s, n) == SettingsDigest(s)));
}

BOOST_AUTO_TEST_CASE(enums_encode_by_fixed_name)
{
	BOOST_CHECK_EQUAL(CanonicalEncoding(tracking::Predictor::RKF45),
	                  "(cfg Predictor choice=RKF45)");
	algorithm::MetaConfig m;
	BOOST_CHECK_EQUAL(CanonicalEncoding(m), "(cfg Meta tracktype=ZeroDim)");
}

// ---- the golden fixture: cross-session digest stability ----

BOOST_AUTO_TEST_CASE(golden_digests_match_committed_fixture)
{
	// Recipes with every field EXPLICITLY set (some structs have uninitialized or
	// precision-dependent defaults; a fixture must never depend on those).
	std::map<std::string, std::string> actual;

	{
		tracking::SteppingConfig c;   // defaults are fully specified for this struct
		actual["stepping_default"] = ConfigDigest(c).Hex();
	}
	{
		tracking::NewtonConfig c;
		c.max_num_newton_iterations = 3;
		c.min_num_newton_iterations = 1;
		actual["newton_3_1"] = ConfigDigest(c).Hex();
	}
	{
		tracking::AdaptiveMultiplePrecisionConfig c;   // defaults; epsilon/Phi/Psi zero-init
		c.epsilon = 4.0; c.Phi = 20000.0; c.Psi = 5000.0;
		actual["amp_explicit"] = ConfigDigest(c).Hex();
	}
	{
		endgame::EndgameConfig c;     // defaults fully specified
		actual["endgame_default"] = ConfigDigest(c).Hex();
	}
	{
		endgame::CauchyConfig c;      // defaults fully specified
		actual["cauchy_default"] = ConfigDigest(c).Hex();
	}
	{
		algorithm::TolerancesConfig c;  // defaults fully specified
		actual["tolerances_default"] = ConfigDigest(c).Hex();
	}
	{
		algorithm::ZeroDimConfig c;
		c.initial_ambient_precision = 30;   // default is DefaultPrecision(): pin explicitly
		actual["zerodim_prec30"] = ConfigDigest(c).Hex();
	}
	{
		tracking::SteppingConfig s;
		tracking::NewtonConfig n;
		endgame::EndgameConfig e;
		actual["composite_stepping_newton_endgame"] = SettingsDigest(s, n, e).Hex();
	}

	auto const fixture_path =
		std::filesystem::path(__FILE__).parent_path() / "data" / "config_digest_fixture.txt";

	std::map<std::string, std::string> expected;
	{
		std::ifstream fixture(fixture_path);
		BOOST_REQUIRE_MESSAGE(fixture.good(), "golden fixture missing: " << fixture_path);
		std::string name, hex;
		while (fixture >> name >> hex)
			expected[name] = hex;
	}

	for (auto const& [name, hex] : actual)
	{
		auto const it = expected.find(name);
		if (it == expected.end())
		{
			BOOST_ERROR("fixture is missing recipe '" << name << "'; its current digest is " << hex);
			continue;
		}
		BOOST_CHECK_MESSAGE(hex == it->second,
			"DIGEST DRIFT for '" << name << "': expected " << it->second << " got " << hex
			<< " -- if intentional, bump b2cfgenc/<n> and the fixture in the same commit");
	}
}

BOOST_AUTO_TEST_SUITE_END()
