//This file is part of Bertini 2.
//
//system_identity_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_identity_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system_identity_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for System content identity (ADR-0042): ContentDigest / Hash / IsSame
semantics, what is and is not identity, and the CROSS-SESSION GOLDEN FIXTURE
(data/system_digest_fixture.txt).  If the golden digests change, the canonical
encoding drifted: that is a digest-breaking event -- bump `b2sysenc/<n>` and
regenerate the fixture in the same commit, never silently.
*/

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"
#include "bertini2/function_tree/canonical.hpp"

using bertini::System;

namespace {

System Parse(std::string const& str)
{
	System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	return sys;
}

std::string const kCircleLine = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";

// Pin the session canonicalization settings for the golden fixture (they are part of the
// encoding header, so unpinned settings would honestly change the digest -- the fixture
// pins one configuration).  RAII-restores the previous settings.
struct PinnedCanonicalization
{
	/// The monomial order in effect before pinning (restored at scope exit).
	bertini::node::MonomialOrder old_order;
	/// The canonicalize-by-default switch before pinning.
	bool old_canon;
	/// The power-fold-by-default switch before pinning.
	bool old_fold;

	PinnedCanonicalization()
		: old_order(bertini::node::CurrentMonomialOrder()),
		  old_canon(bertini::node::CanonicalizeByDefault()),
		  old_fold(bertini::node::PowerFoldByDefault())
	{
		bertini::node::SetMonomialOrder(bertini::node::MonomialOrder::GrevLex);
		bertini::node::SetCanonicalizeByDefault(true);
		bertini::node::SetPowerFoldByDefault(true);
	}
	~PinnedCanonicalization()
	{
		bertini::node::SetMonomialOrder(old_order);
		bertini::node::SetCanonicalizeByDefault(old_canon);
		bertini::node::SetPowerFoldByDefault(old_fold);
	}
};

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(system_identity)

// ---- IsSame / Hash semantics ----

BOOST_AUTO_TEST_CASE(is_same_is_reflexive_and_symmetric)
{
	auto a = Parse(kCircleLine);
	auto b = Parse(kCircleLine);   // independently built, equal content
	BOOST_CHECK(a.IsSame(a));
	BOOST_CHECK(a.IsSame(b));
	BOOST_CHECK(b.IsSame(a));
	BOOST_CHECK_EQUAL(a.Hash(), b.Hash());
	BOOST_CHECK_EQUAL(a.ContentDigest().Hex(), b.ContentDigest().Hex());
}

BOOST_AUTO_TEST_CASE(different_functions_differ)
{
	auto a = Parse(kCircleLine);
	auto b = Parse("function f,g; variable_group x,y; f = x^2+y^2-2; g = x-y;");
	BOOST_CHECK(!a.IsSame(b));
}

BOOST_AUTO_TEST_CASE(different_variable_grouping_differs)
{
	auto grouped_together = Parse("function f; variable_group x,y; f = x*y-1;");
	auto grouped_apart = Parse("function f; variable_group x; variable_group y; f = x*y-1;");
	BOOST_CHECK(!grouped_together.IsSame(grouped_apart));
}

BOOST_AUTO_TEST_CASE(variable_names_are_identity)
{
	auto x_system = Parse("function f; variable_group x; f = x^2-1;");
	auto z_system = Parse("function f; variable_group z; f = z^2-1;");
	BOOST_CHECK(!x_system.IsSame(z_system));
}

// ---- transient state is NOT identity ----

BOOST_AUTO_TEST_CASE(precision_and_differentiation_do_not_change_digest)
{
	auto sys = Parse(kCircleLine);
	auto const before = sys.ContentDigest();

	sys.precision(50);
	BOOST_CHECK(before == sys.ContentDigest());

	sys.Differentiate();
	BOOST_CHECK(before == sys.ContentDigest());
}

// ---- structural operations ARE identity ----

BOOST_AUTO_TEST_CASE(homogenize_changes_digest)
{
	auto sys = Parse(kCircleLine);
	auto const before = sys.ContentDigest();
	sys.Homogenize();
	BOOST_CHECK(!(before == sys.ContentDigest()));
}

BOOST_AUTO_TEST_CASE(autopatch_changes_digest_and_random_patches_differ)
{
	auto a = Parse(kCircleLine);
	auto b = Parse(kCircleLine);
	a.Homogenize(); b.Homogenize();
	BOOST_CHECK(a.IsSame(b));

	auto const unpatched = a.ContentDigest();
	a.AutoPatch(); b.AutoPatch();
	BOOST_CHECK(!(unpatched == a.ContentDigest()));  // patching is identity
	BOOST_CHECK(!a.IsSame(b));                        // two random patches differ
}

BOOST_AUTO_TEST_CASE(random_start_systems_differ_and_gamma_is_identity)
{
	auto sys = Parse(kCircleLine);
	sys.Homogenize();
	sys.AutoPatch();

	// two total-degree starts draw different random coefficients: different identity
	bertini::start_system::TotalDegreeBinomial start_a(sys), start_b(sys);
	BOOST_CHECK(!start_a.IsSame(start_b));

	// two homotopies over the SAME start differ only by gamma -- and gamma is identity
	auto t = bertini::node::Variable::Make("t");
	auto gamma_a = bertini::node::Rational::Make(3, 7, 1, 11);
	auto gamma_b = bertini::node::Rational::Make(5, 13, 2, 9);
	auto homotopy_a = (1-t)*sys + gamma_a*t*start_a;
	auto homotopy_b = (1-t)*sys + gamma_b*t*start_a;
	homotopy_a.AddPathVariable(t);
	homotopy_b.AddPathVariable(t);
	BOOST_CHECK(!homotopy_a.IsSame(homotopy_b));

	// same gamma, same start: same identity
	auto homotopy_c = (1-t)*sys + gamma_a*t*start_a;
	homotopy_c.AddPathVariable(t);
	BOOST_CHECK(homotopy_a.IsSame(homotopy_c));
}

BOOST_AUTO_TEST_CASE(serialization_round_trip_preserves_digest)
{
	auto original = Parse(kCircleLine);
	original.Homogenize();
	original.AutoPatch();

	std::stringstream archive_stream;
	{
		boost::archive::text_oarchive oa(archive_stream);
		oa << original;
	}
	System loaded;
	{
		boost::archive::text_iarchive ia(archive_stream);
		ia >> loaded;
	}
	BOOST_CHECK_EQUAL(original.ContentDigest().Hex(), loaded.ContentDigest().Hex());
	BOOST_CHECK(original.IsSame(loaded));
}

// ---- Seal(): hashcons-on-freeze ----

BOOST_AUTO_TEST_CASE(copy_is_same_as_original_including_prehom_snapshot)
{
	auto original = Parse(kCircleLine);
	original.Homogenize();   // populates the pre-homogenization snapshot
	System const copy(original);
	BOOST_CHECK(copy.IsSame(original));   // regression: copy ctor must carry prehom functions
}

BOOST_AUTO_TEST_CASE(seal_memoizes_and_is_idempotent)
{
	auto sys = Parse(kCircleLine);
	auto const fresh = sys.ContentDigest();
	BOOST_CHECK(!sys.IsSealed());

	sys.Seal();
	BOOST_CHECK(sys.IsSealed());
	BOOST_CHECK(fresh == sys.ContentDigest());   // memo agrees with fresh compute
	sys.Seal();                                  // idempotent
	BOOST_CHECK(fresh == sys.ContentDigest());
}

BOOST_AUTO_TEST_CASE(every_structural_mutator_throws_when_sealed)
{
	auto sys = Parse(kCircleLine);
	sys.Seal();

	auto x = bertini::node::Variable::Make("x");
	auto w = bertini::node::Variable::Make("w");
	auto other = Parse(kCircleLine);

	BOOST_CHECK_THROW(sys.AddFunction(x + 1), std::logic_error);
	BOOST_CHECK_THROW(sys.AddFunctions({x + 1}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddVariableGroup({w}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddHomVariableGroup({w}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddUngroupedVariable(w), std::logic_error);
	BOOST_CHECK_THROW(sys.AddUngroupedVariables({w}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddImplicitParameter(w), std::logic_error);
	BOOST_CHECK_THROW(sys.AddImplicitParameters({w}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddPathVariable(w), std::logic_error);
	BOOST_CHECK_THROW(sys.Homogenize(), std::logic_error);
	BOOST_CHECK_THROW(sys.AutoPatch(), std::logic_error);
	BOOST_CHECK_THROW(sys.CopyPatches(other), std::logic_error);
	BOOST_CHECK_THROW(sys.CopyVariableStructure(other), std::logic_error);
	BOOST_CHECK_THROW(sys.ClearFunctions(), std::logic_error);
	BOOST_CHECK_THROW(sys.ClearBlocks(), std::logic_error);
	BOOST_CHECK_THROW(sys.ClearVariables(), std::logic_error);
	BOOST_CHECK_THROW(sys.SimplifyFunctions(), std::logic_error);
	BOOST_CHECK_THROW(sys.Simplify(), std::logic_error);
	BOOST_CHECK_THROW(sys.ReorderFunctionsByDegreeDecreasing(), std::logic_error);
	BOOST_CHECK_THROW(sys.ReorderFunctionsByDegreeIncreasing(), std::logic_error);
	BOOST_CHECK_THROW(sys += other, std::logic_error);
	BOOST_CHECK_THROW(sys *= (x + 1), std::logic_error);
	BOOST_CHECK_THROW(sys.SetVariableGroups({{w}}), std::logic_error);
	BOOST_CHECK_THROW(sys.AddBlock(bertini::blocks::PolynomialBlock{}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(transient_operations_work_on_a_sealed_system)
{
	auto sys = Parse(kCircleLine);
	sys.Seal();
	auto const digest = sys.ContentDigest();

	sys.precision(50);      // transient
	sys.Differentiate();    // derived cache

	// evaluation on a sealed system
	bertini::Vec<bertini::complex_dbl> values(2), point(2);
	point << bertini::complex_dbl(1.0, 0.0), bertini::complex_dbl(0.0, 1.0);
	sys.EvalInPlace(values, point);

	BOOST_CHECK(digest == sys.ContentDigest());
	BOOST_CHECK(sys.IsSealed());
}

BOOST_AUTO_TEST_CASE(copy_of_sealed_system_is_unsealed_and_mutable)
{
	auto sys = Parse(kCircleLine);
	sys.Seal();

	System copy(sys);
	BOOST_CHECK(!copy.IsSealed());
	BOOST_CHECK(copy.IsSame(sys));

	auto w = bertini::node::Variable::Make("w");
	BOOST_CHECK_NO_THROW(copy.AddUngroupedVariable(w));   // the copy-on-write escape hatch
	BOOST_CHECK(!copy.IsSame(sys));

	System assigned;
	assigned = sys;
	BOOST_CHECK(!assigned.IsSealed());
}

BOOST_AUTO_TEST_CASE(seal_flag_round_trips_serialization_and_digest_rememoizes)
{
	auto original = Parse(kCircleLine);
	original.Seal();

	std::stringstream archive_stream;
	{
		boost::archive::text_oarchive oa(archive_stream);
		oa << original;
	}
	System loaded;
	{
		boost::archive::text_iarchive ia(archive_stream);
		ia >> loaded;
	}
	BOOST_CHECK(loaded.IsSealed());
	BOOST_CHECK(loaded.IsSame(original));   // digest recomputed on demand, then memoized
	auto w = bertini::node::Variable::Make("w");
	BOOST_CHECK_THROW(loaded.AddUngroupedVariable(w), std::logic_error);
}

// ---- InternSystem: the weak intern table ----

BOOST_AUTO_TEST_CASE(intern_hit_returns_the_live_representative)
{
	auto first = std::make_shared<System>(Parse(kCircleLine));
	auto second = std::make_shared<System>(Parse(kCircleLine));   // equal content, distinct object
	BOOST_REQUIRE(first != second);

	auto rep_a = bertini::InternSystem(first);
	BOOST_CHECK(rep_a == first);            // miss: candidate becomes representative
	BOOST_CHECK(first->IsSealed());         // interning seals

	auto rep_b = bertini::InternSystem(second);
	BOOST_CHECK(rep_b == rep_a);            // hit: the SAME shared object comes back
}

BOOST_AUTO_TEST_CASE(intern_distinguishes_different_content)
{
	auto a = std::make_shared<System>(Parse("function f; variable_group x; f = x^5-1;"));
	auto b = std::make_shared<System>(Parse("function f; variable_group x; f = x^5-2;"));
	BOOST_CHECK(bertini::InternSystem(a) != bertini::InternSystem(b));
}

BOOST_AUTO_TEST_CASE(intern_table_is_weak_and_self_cleans)
{
	std::string const recipe = "function f; variable_group x; f = x^7 - 42;";
	{
		auto ephemeral = std::make_shared<System>(Parse(recipe));
		auto rep = bertini::InternSystem(ephemeral);
		BOOST_CHECK(rep == ephemeral);
	}   // last owner drops; the table's weak_ptr is now expired

	// a fresh equal system becomes a NEW representative (no dangling resurrection)
	auto fresh = std::make_shared<System>(Parse(recipe));
	auto rep = bertini::InternSystem(fresh);
	BOOST_CHECK(rep == fresh);
}

// ---- the golden fixture: cross-session digest stability ----

BOOST_AUTO_TEST_CASE(golden_digests_match_committed_fixture)
{
	PinnedCanonicalization const pin;

	// deterministic recipes only (no randomness), so the digests are reproducible forever
	std::map<std::string, System> recipes;
	recipes.emplace("single_var", Parse("function f; variable_group x; f = x+1;"));
	recipes.emplace("circle_line", Parse(kCircleLine));
	recipes.emplace("rational_coeffs", Parse("function f; variable_group x,y; f = (1/3)*x^3 - y + 2;"));
	recipes.emplace("transcendental", Parse("function f; variable_group x; f = sin(x) + 3*exp(x);"));
	{
		auto homogenized = Parse(kCircleLine);
		homogenized.Homogenize();   // deterministic: no random coefficients involved
		recipes.emplace("circle_line_homogenized", std::move(homogenized));
	}

	auto const fixture_path =
		std::filesystem::path(__FILE__).parent_path() / "data" / "system_digest_fixture.txt";

	std::map<std::string, std::string> expected;
	{
		std::ifstream fixture(fixture_path);
		BOOST_REQUIRE_MESSAGE(fixture.good(), "golden fixture missing: " << fixture_path);
		std::string name, hex;
		while (fixture >> name >> hex)
			expected[name] = hex;
	}

	for (auto const& [name, sys] : recipes)
	{
		auto const it = expected.find(name);
		if (it == expected.end())
		{
			BOOST_ERROR("fixture is missing recipe '" << name << "'; its current digest is "
				<< sys.ContentDigest().Hex());
			continue;
		}
		BOOST_CHECK_MESSAGE(sys.ContentDigest().Hex() == it->second,
			"DIGEST DRIFT for '" << name << "': expected " << it->second
			<< " got " << sys.ContentDigest().Hex()
			<< " -- if intentional, bump b2sysenc/<n> and the fixture in the same commit");
	}
}

BOOST_AUTO_TEST_SUITE_END()
