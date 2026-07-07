//This file is part of Bertini 2.
//
//reintern_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//reintern_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with reintern_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for node::Reintern (ADR-0042): deserialized DAGs rebuilt through the
interning factories unify pointer-wise with live equivalents, preserve the archive's
internal sharing, and never change content (canonical encoding is bit-identical).
*/

#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical_encoding.hpp"
#include "bertini2/function_tree/reintern.hpp"

using namespace bertini::node;
using Nd = std::shared_ptr<Node>;

namespace {

// serialize -> deserialize a tree, yielding a raw (un-interned) clone
Nd RoundTrip(Nd const& tree)
{
	std::stringstream archive_stream;
	{
		boost::archive::text_oarchive oa(archive_stream);
		oa << tree;
	}
	Nd loaded;
	{
		boost::archive::text_iarchive ia(archive_stream);
		ia >> loaded;
	}
	return loaded;
}

Nd BuildMixedTree()
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	return pow(x, 3) + Rational::Make(1, 3, 0, 1) * y - sin(x) * Pi() + exp(y) / 2;
}

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(reintern)

BOOST_AUTO_TEST_CASE(loaded_tree_does_not_share_until_reinterned)
{
	auto live = BuildMixedTree();
	auto loaded = RoundTrip(live);

	BOOST_CHECK(loaded != live);                    // deserialization bypasses Intern
	BOOST_CHECK(Reintern(loaded) == live);          // reinterning unifies pointer-wise
}

BOOST_AUTO_TEST_CASE(reintern_never_changes_content)
{
	auto live = BuildMixedTree();
	auto loaded = RoundTrip(live);
	auto reinterned = Reintern(loaded);
	BOOST_CHECK_EQUAL(CanonicalEncoding(reinterned), CanonicalEncoding(live));
}

BOOST_AUTO_TEST_CASE(variables_unify_by_name)
{
	auto loaded = RoundTrip(Variable::Make("x"));
	BOOST_CHECK(Reintern(loaded) == Variable::Make("x"));
}

BOOST_AUTO_TEST_CASE(two_archives_unify_with_each_other)
{
	auto live = BuildMixedTree();
	auto first = Reintern(RoundTrip(live));
	auto second = Reintern(RoundTrip(live));
	BOOST_CHECK(first == second);                   // both collapse to one interned object
}

BOOST_AUTO_TEST_CASE(memo_preserves_archive_sharing_across_roots)
{
	auto x = Variable::Make("x");
	auto shared_term = pow(x, 2) + 1;
	Nd f = shared_term + x;
	Nd g = shared_term - 5;

	// archive both roots together: boost object tracking preserves their sharing
	std::stringstream archive_stream;
	{
		boost::archive::text_oarchive oa(archive_stream);
		oa << f << g;
	}
	Nd loaded_f, loaded_g;
	{
		boost::archive::text_iarchive ia(archive_stream);
		ia >> loaded_f >> loaded_g;
	}

	ReinternMemo memo;
	auto ref = Reintern(loaded_f, memo);
	auto reg = Reintern(loaded_g, memo);

	BOOST_CHECK(ref == f);
	BOOST_CHECK(reg == g);
	// and the shared subtree is rebuilt exactly once (memo hit): both rebuilt roots hold
	// the SAME interned shared_term the live trees do -- checked via pointer unification
	// of the whole roots above.
	BOOST_CHECK_EQUAL(memo.count(loaded_f.get()), 1u);
	BOOST_CHECK_EQUAL(memo.count(loaded_g.get()), 1u);
}

BOOST_AUTO_TEST_CASE(every_kind_survives_reintern)
{
	auto x = Variable::Make("x");
	std::vector<Nd> const kinds = {
		Variable::Make("v"), Integer::Make(-3), Rational::Make(2, 7, 1, 5),
		Complex::Make("1.25", "-0.5"), Pi(), E(),
		x + 1, x * 2, -x, pow(x, x), pow(x, 4), sqrt(x), exp(x), log(x),
		sin(x), asin(x), cos(x), acos(x), tan(x), atan(x),
		Named(x + 2, "fun")
	};
	for (auto const& k : kinds)
	{
		auto loaded = RoundTrip(k);
		auto reinterned = Reintern(loaded);
		BOOST_CHECK_MESSAGE(reinterned == k,
			"kind failed to unify: " << CanonicalEncoding(k));
		BOOST_CHECK_EQUAL(CanonicalEncoding(reinterned), CanonicalEncoding(k));
	}
}

BOOST_AUTO_TEST_SUITE_END()
