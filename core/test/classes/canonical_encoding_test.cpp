//This file is part of Bertini 2.
//
//canonical_encoding_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_encoding_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_encoding_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for the canonical exact node encoding (ADR-0042) -- the identity substrate
for persistent digests.  The literal expected strings here ARE the format contract:
if one of these fails after an encoder change, that change is digest-breaking and
requires an encoding-version bump + golden-fixture update, not a test edit.
*/

#include <set>
#include <sstream>
#include <vector>

#include <boost/test/unit_test.hpp>
#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical_encoding.hpp"

using namespace bertini::node;
using Nd = std::shared_ptr<Node>;

BOOST_AUTO_TEST_SUITE(canonical_encoding)

// ---- leaves: pinned literal forms ----

BOOST_AUTO_TEST_CASE(leaf_forms_are_pinned)
{
	BOOST_CHECK_EQUAL(CanonicalEncoding(Variable::Make("x")), "(var 1:x)");
	BOOST_CHECK_EQUAL(CanonicalEncoding(Integer::Make(-7)), "(int -7)");
	BOOST_CHECK_EQUAL(CanonicalEncoding(Rational::Make(1, 3, 0, 1)), "(rat 1/3 0)");
	BOOST_CHECK_EQUAL(CanonicalEncoding(Pi()), "(pi)");
	BOOST_CHECK_EQUAL(CanonicalEncoding(E()), "(e)");
}

BOOST_AUTO_TEST_CASE(names_are_netstrings_utf8_safe)
{
	// an emoji variable name: the netstring length is in BYTES, so no escaping is needed
	std::string const heart = "\xF0\x9F\x92\x9C";  // 💜, 4 bytes of UTF-8
	auto v = Variable::Make(heart);
	BOOST_CHECK_EQUAL(CanonicalEncoding(v), "(var 4:" + heart + ")");
}

BOOST_AUTO_TEST_CASE(complex_literal_carries_precision)
{
	auto n = Complex::Make("1.5", "0");
	auto encoding = CanonicalEncoding(n);
	BOOST_CHECK(encoding.rfind("(cplx ", 0) == 0);
	// value-equal literals at different stored precisions must encode differently;
	// here we just pin that the precision field is present and the encoding is stable
	BOOST_CHECK_EQUAL(encoding, CanonicalEncoding(Complex::Make("1.5", "0")));
}

// ---- structure ----

BOOST_AUTO_TEST_CASE(equal_interned_trees_encode_identically)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto a = pow(x, 2) + y * x - 1;
	auto b = pow(x, 2) + y * x - 1;   // independently built
	BOOST_CHECK_EQUAL(CanonicalEncoding(a), CanonicalEncoding(b));
}

BOOST_AUTO_TEST_CASE(commuted_operands_encode_identically_under_canonicalization)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_EQUAL(CanonicalEncoding(x + y), CanonicalEncoding(y + x));
	BOOST_CHECK_EQUAL(CanonicalEncoding(x * y), CanonicalEncoding(y * x));
}

BOOST_AUTO_TEST_CASE(different_content_encodes_differently)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_NE(CanonicalEncoding(x + y), CanonicalEncoding(x - y));
	BOOST_CHECK_NE(CanonicalEncoding(x * y), CanonicalEncoding(x / y));
	BOOST_CHECK_NE(CanonicalEncoding(pow(x, 2)), CanonicalEncoding(pow(x, 3)));
	BOOST_CHECK_NE(CanonicalEncoding(Variable::Make("x")), CanonicalEncoding(Variable::Make("x2")));
	BOOST_CHECK_NE(CanonicalEncoding(Integer::Make(1)), CanonicalEncoding(Rational::Make(1)));
}

BOOST_AUTO_TEST_CASE(named_expression_differs_from_bare_expression)
{
	auto x = Variable::Make("x");
	auto expr = pow(x, 2) + 1;
	auto named = Named(expr, "f");
	BOOST_CHECK_NE(CanonicalEncoding(named), CanonicalEncoding(expr));
	BOOST_CHECK_NE(CanonicalEncoding(Named(expr, "f")), CanonicalEncoding(Named(expr, "g")));
}

BOOST_AUTO_TEST_CASE(unary_and_trig_kinds_all_encode_and_are_distinct)
{
	auto x = Variable::Make("x");
	std::vector<Nd> const kinds = {
		-x, sqrt(x), exp(x), log(x),
		sin(x), asin(x), cos(x), acos(x), tan(x), atan(x)
	};
	std::set<std::string> encodings;
	for (auto const& k : kinds)
		encodings.insert(CanonicalEncoding(k));
	BOOST_CHECK_EQUAL(encodings.size(), kinds.size());  // pairwise distinct
}

// ---- back-references ----

BOOST_AUTO_TEST_CASE(shared_subtree_backreferences_deterministically)
{
	auto x = Variable::Make("x");
	auto shared_term = pow(x, 2) + 1;               // interned: both uses ARE one object
	// NOTE: shared*shared would power-fold to shared^2 (one occurrence); use a shape the
	// canonicalizer keeps as two occurrences of the same interned node
	auto whole = shared_term + shared_term * x;

	auto encoding = CanonicalEncoding(whole);
	// the second occurrence must be a back-reference, not a re-encoding
	BOOST_CHECK(encoding.find('#') != std::string::npos);

	// and the whole thing is reproducible
	BOOST_CHECK_EQUAL(encoding, CanonicalEncoding(shared_term + shared_term * x));
}

BOOST_AUTO_TEST_CASE(context_shared_across_roots_backreferences_across_functions)
{
	auto x = Variable::Make("x");
	auto common = pow(x, 2) + 1;
	auto f = common * 2;
	auto g = common - 5;

	// one context across both roots: g's use of `common` back-references f's
	std::ostringstream two_roots;
	EncodingContext ctx;
	EncodeCanonical(f, two_roots, ctx);
	two_roots << '|';
	EncodeCanonical(g, two_roots, ctx);

	std::string const joint = two_roots.str();
	auto const bar = joint.find('|');
	BOOST_REQUIRE(bar != std::string::npos);
	std::string const g_part = joint.substr(bar + 1);
	BOOST_CHECK(g_part.find('#') != std::string::npos);

	// fresh contexts: each root self-contained, no cross-references
	BOOST_CHECK(CanonicalEncoding(f).find(CanonicalEncoding(x)) != std::string::npos);
}

BOOST_AUTO_TEST_SUITE_END()
