//This file is part of Bertini 2.
//
//canonical_decoding_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_decoding_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_decoding_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical_decoding_test.cpp

\brief The canonical-encoding READER (ADR-0042): every node kind and every system section
round-trips, and the proof is content identity -- the decoded object encodes to the same
text and carries the same digest.  Malformed, foreign-version and foreign-settings texts
are refused loudly.
*/

#include <boost/test/unit_test.hpp>

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical.hpp"
#include "bertini2/function_tree/canonical_encoding.hpp"
#include "bertini2/function_tree/canonical_decoding.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

using namespace bertini;
using namespace bertini::node;
using Nd = std::shared_ptr<Node>;

namespace {

System Parse(std::string const& str)
{
	System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	return sys;
}

std::string const kCircleLine = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";

// Pin the session canonicalization settings, as the golden-fixture test does; restored at exit.
struct PinnedCanonicalization
{
	/// The monomial order in effect before pinning (restored at scope exit).
	MonomialOrder old_order;
	/// The canonicalize-by-default switch before pinning.
	bool old_canon;
	/// The power-fold-by-default switch before pinning.
	bool old_fold;

	PinnedCanonicalization()
		: old_order(CurrentMonomialOrder()), old_canon(CanonicalizeByDefault()), old_fold(PowerFoldByDefault())
	{
		SetMonomialOrder(MonomialOrder::GrevLex);
		SetCanonicalizeByDefault(true);
		SetPowerFoldByDefault(true);
	}
	~PinnedCanonicalization()
	{
		SetMonomialOrder(old_order);
		SetCanonicalizeByDefault(old_canon);
		SetPowerFoldByDefault(old_fold);
	}
};

// A complex literal at a chosen precision, digits taken exactly.
complex_mp ComplexAt(unsigned precision, std::string const& re, std::string const& im)
{
	real_mp const r(re, precision);
	real_mp const i(im, precision);
	complex_mp z;
	z.precision(precision);
	z.real(r);
	z.imag(i);
	return z;
}

// The round trip for a node: decode(encode(n)) encodes identically, and -- the DAG being
// interned -- IS the same object.
void CheckNodeRoundTrip(Nd const& n)
{
	auto const text = CanonicalEncoding(n);
	auto decoded = DecodeCanonicalTree(text);
	BOOST_CHECK_EQUAL(CanonicalEncoding(decoded), text);
	BOOST_CHECK(decoded == n);
}

// The round trip for a system: same text, same digest, same values.
void CheckSystemRoundTrip(System const& sys)
{
	auto const text = sys.CanonicalEncodingText();
	auto back = System::FromCanonicalEncoding(text);
	BOOST_CHECK_EQUAL(back.CanonicalEncodingText(), text);
	BOOST_CHECK(back.ContentDigest() == sys.ContentDigest());
	BOOST_CHECK(back.IsSame(sys));

	BOOST_REQUIRE_EQUAL(back.NumVariables(), sys.NumVariables());
	BOOST_REQUIRE_EQUAL(back.NumTotalFunctions(), sys.NumTotalFunctions());
	Vec<complex_dbl> point(sys.NumVariables());
	for (Eigen::Index ii = 0; ii < point.size(); ++ii)
		point(ii) = complex_dbl(0.3 + 0.1 * double(ii), 0.7 - 0.2 * double(ii));
	Vec<complex_dbl> expected, actual;
	if (sys.HavePathVariable())
	{
		expected = sys.Eval(point, complex_dbl(0.6, 0.1));
		actual = back.Eval(point, complex_dbl(0.6, 0.1));
	}
	else
	{
		expected = sys.Eval(point);
		actual = back.Eval(point);
	}
	BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
	BOOST_CHECK_SMALL((actual - expected).norm(), 1e-12 * (1 + expected.norm()));
}

} // unnamed namespace


BOOST_AUTO_TEST_SUITE(canonical_decoding)

// ---- nodes: every kind the encoder writes, the reader reads ----

BOOST_AUTO_TEST_CASE(leaves_round_trip)
{
	CheckNodeRoundTrip(Variable::Make("x"));
	CheckNodeRoundTrip(Integer::Make(-7));
	CheckNodeRoundTrip(Integer::Make(mpz_int("123456789012345678901234567890")));
	CheckNodeRoundTrip(Rational::Make(1, 3, 0, 1));
	CheckNodeRoundTrip(Rational::Make(-3, 7, 5, 11));
	CheckNodeRoundTrip(Pi());
	CheckNodeRoundTrip(E());
}

BOOST_AUTO_TEST_CASE(names_with_any_bytes_round_trip)
{
	std::string const heart = "\xF0\x9F\x92\x9C";   // 💜, 4 bytes: the netstring length is in bytes
	CheckNodeRoundTrip(Variable::Make(heart));
	CheckNodeRoundTrip(Variable::Make("x_" + heart + "1"));
	CheckNodeRoundTrip(NamedExpression::Make(Variable::Make("x") + 1, "f" + heart));
}

BOOST_AUTO_TEST_CASE(complex_literals_keep_their_stored_precision_and_digits)
{
	auto at50 = Complex::Make(ComplexAt(50, "1.23456789012345678901234567890123456789012345678901", "-2.5"));
	auto at20 = Complex::Make(ComplexAt(20, "0.1", "0.2"));
	CheckNodeRoundTrip(at50);
	CheckNodeRoundTrip(at20);

	auto decoded = std::dynamic_pointer_cast<Complex>(DecodeCanonicalTree(CanonicalEncoding(at50)));
	BOOST_REQUIRE(decoded);
	BOOST_CHECK_EQUAL(decoded->GetValue().precision(), 50u);
	BOOST_CHECK(decoded->GetValue() == at50->GetValue());

	// a product mixing precisions keeps each literal's own precision
	auto x = Variable::Make("x");
	CheckNodeRoundTrip(at50 * x + at20 * pow(x, 2));
}

BOOST_AUTO_TEST_CASE(operators_round_trip)
{
	auto x = Variable::Make("x"), y = Variable::Make("y");
	CheckNodeRoundTrip(x + y - 3);                             // signs in a sum
	CheckNodeRoundTrip(x * y / 7);                             // multiply and divide
	CheckNodeRoundTrip(pow(x, 5));                             // integer power
	CheckNodeRoundTrip(pow(x, Rational::Make(1, 2, 0, 1)));    // a non-integer power
	CheckNodeRoundTrip(pow(x + 1, y));                         // a symbolic power
	CheckNodeRoundTrip(NamedExpression::Make(pow(x, 2) + pow(y, 2) - 1, "f"));
}

BOOST_AUTO_TEST_CASE(unary_kinds_round_trip)
{
	auto x = Variable::Make("x");
	std::vector<Nd> const kinds = {
		-x, sqrt(x), exp(x), log(x),
		sin(x), asin(x), cos(x), acos(x), tan(x), atan(x)
	};
	for (auto const& k : kinds)
		CheckNodeRoundTrip(k);
	CheckNodeRoundTrip(sin(exp(x) - 1) * atan(pow(x, 3)));
}

BOOST_AUTO_TEST_CASE(differentials_round_trip)
{
	// differentials are not hash-consed, so the check is text equality, not object identity
	auto x = Variable::Make("x");
	auto dx = Differential::Make(x, "x");
	for (auto const& n : {Nd(dx), Nd(2 * x * dx)})
	{
		auto const text = CanonicalEncoding(n);
		BOOST_CHECK_EQUAL(CanonicalEncoding(DecodeCanonicalTree(text)), text);
	}
}

BOOST_AUTO_TEST_CASE(back_references_resolve_to_the_shared_node)
{
	auto x = Variable::Make("x");
	auto shared_term = pow(x, 2) + 1;
	auto whole = shared_term + shared_term * x;
	auto const text = CanonicalEncoding(whole);
	BOOST_REQUIRE(text.find('#') != std::string::npos);   // the shape does back-reference
	CheckNodeRoundTrip(whole);
}

BOOST_AUTO_TEST_CASE(a_shared_context_resolves_back_references_across_roots)
{
	auto x = Variable::Make("x");
	auto shared_term = pow(x, 3) - x;
	auto f = shared_term + 1;
	auto g = shared_term * 2;

	std::ostringstream out;
	EncodingContext ectx;
	EncodeCanonical(f, out, ectx);
	out << ' ';
	EncodeCanonical(g, out, ectx);            // g's copy of the shared term is a back-reference
	auto const text = out.str();

	DecodingCursor cursor(text);
	DecodingContext dctx;
	auto f_back = DecodeCanonical(cursor, dctx);
	cursor.Expect(' ');
	auto g_back = DecodeCanonical(cursor, dctx);
	BOOST_CHECK(cursor.AtEnd());
	BOOST_CHECK(f_back == f);
	BOOST_CHECK(g_back == g);
}

// ---- systems: every section, every block kind ----

BOOST_AUTO_TEST_CASE(plain_polynomial_systems_round_trip)
{
	PinnedCanonicalization const pin;
	CheckSystemRoundTrip(Parse("function f; variable_group x; f = x+1;"));
	CheckSystemRoundTrip(Parse(kCircleLine));
	CheckSystemRoundTrip(Parse("function f; variable_group x,y; f = (1/3)*x^3 - y + 2;"));
	CheckSystemRoundTrip(Parse("function f; variable_group x; f = sin(x) + 3*exp(x);"));
	CheckSystemRoundTrip(Parse("function f,g; variable_group x; variable_group y; f = x*y - 1; g = x - 2*y;"));
	CheckSystemRoundTrip(Parse("function f; hom_variable_group x,y,z; f = x^2 + y*z;"));
}

BOOST_AUTO_TEST_CASE(homogenized_and_patched_systems_round_trip)
{
	PinnedCanonicalization const pin;
	auto homogenized = Parse(kCircleLine);
	homogenized.Homogenize();
	CheckSystemRoundTrip(homogenized);       // homogenizing variables + the pre-homogenization snapshot

	auto patched = Parse(kCircleLine);
	patched.Homogenize();
	patched.AutoPatch();                     // random patch coefficients, exact in the encoding
	CheckSystemRoundTrip(patched);
	auto back = System::FromCanonicalEncoding(patched.CanonicalEncodingText());
	BOOST_CHECK(back.IsPatched());
}

BOOST_AUTO_TEST_CASE(a_homotopy_with_a_blend_block_round_trips)
{
	PinnedCanonicalization const pin;
	auto target = Parse(kCircleLine);
	target.Homogenize();
	target.AutoPatch();
	bertini::start_system::TotalDegreeBinomial start(target);
	auto t = Variable::Make("t");
	auto homotopy = (1-t)*target + Rational::Make(3, 7, 1, 11)*t*start;
	homotopy.AddPathVariable(t);
	BOOST_REQUIRE(homotopy.HavePathVariable());
	CheckSystemRoundTrip(homotopy);

	auto back = System::FromCanonicalEncoding(homotopy.CanonicalEncodingText());
	BOOST_CHECK(back.GetPathVariable() == homotopy.GetPathVariable());   // the interned variable
}

BOOST_AUTO_TEST_CASE(a_randomized_system_round_trips)
{
	PinnedCanonicalization const pin;
	auto over = Parse("function f,g,h; variable_group x,y; f = x^2+y^2-1; g = x-y; h = x*y - 1/4;");
	auto randomized = over.Randomize();     // a randomization block with an operand system inside
	CheckSystemRoundTrip(randomized);

	auto homogenized = over.Randomize();
	homogenized.Homogenize();               // the block records its homogenizing variables
	CheckSystemRoundTrip(homogenized);
}

BOOST_AUTO_TEST_CASE(linear_forms_and_products_of_linears_blocks_round_trip)
{
	PinnedCanonicalization const pin;
	auto x = Variable::Make("x"), y = Variable::Make("y");

	Mat<complex_mp> forms(2, 3);
	forms << ComplexAt(30, "1.5", "-0.25"), ComplexAt(50, "0.1", "0.2"), ComplexAt(30, "-3", "0"),
	         ComplexAt(30, "2", "1"), ComplexAt(30, "0", "-1"), ComplexAt(40, "0.333333333333333333333333333333333333", "0");
	System sliced;
	sliced.AddVariableGroup(VariableGroup{x, y});
	sliced.AddFunction(pow(x, 2) + pow(y, 2) - 1);
	sliced.AddBlock(blocks::LinearFormsBlock(2, forms));
	CheckSystemRoundTrip(sliced);

	Mat<complex_mp> factor_a(2, 3), factor_b(1, 3);
	factor_a << ComplexAt(30, "1", "0"), ComplexAt(30, "2", "0"), ComplexAt(30, "-1", "0.5"),
	            ComplexAt(30, "0", "1"), ComplexAt(30, "1", "1"), ComplexAt(30, "0.75", "0");
	factor_b << ComplexAt(30, "3", "0"), ComplexAt(30, "-1", "0"), ComplexAt(30, "0.5", "-0.5");
	System products;
	products.AddVariableGroup(VariableGroup{x, y});
	products.AddBlock(blocks::ProductsOfLinearsBlock(2, std::vector<Mat<complex_mp>>{factor_a, factor_b}));
	CheckSystemRoundTrip(products);
}

BOOST_AUTO_TEST_CASE(decoded_systems_share_the_live_interned_nodes)
{
	PinnedCanonicalization const pin;
	auto sys = Parse(kCircleLine);
	auto back = System::FromCanonicalEncoding(sys.CanonicalEncodingText());
	BOOST_REQUIRE_EQUAL(back.NumVariables(), 2u);
	BOOST_CHECK(back.Variables()[0] == sys.Variables()[0]);   // the same Variable object, by name
	BOOST_CHECK(back.Variables()[1] == sys.Variables()[1]);
}

// ---- refusals: a wrong text never becomes a wrong system ----

BOOST_AUTO_TEST_CASE(a_foreign_version_is_refused)
{
	PinnedCanonicalization const pin;
	auto text = Parse(kCircleLine).CanonicalEncodingText();
	BOOST_REQUIRE_EQUAL(text.rfind(SystemEncodingVersion, 0), 0u);
	text.replace(0, std::string(SystemEncodingVersion).size(), "b2sysenc/999");
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(text), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(other_canonicalization_settings_are_refused)
{
	PinnedCanonicalization const pin;
	auto const text = Parse(kCircleLine).CanonicalEncodingText();
	SetCanonicalizeByDefault(false);          // the pin restores this at scope exit
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(text), std::runtime_error);
	SetCanonicalizeByDefault(true);
	SetMonomialOrder(MonomialOrder::Lex);
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(text), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(damaged_text_is_refused)
{
	PinnedCanonicalization const pin;
	auto const text = Parse(kCircleLine).CanonicalEncodingText();
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(text.substr(0, text.size() / 2)), std::runtime_error);   // truncated
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(text + "extra"), std::runtime_error);                   // trailing bytes
	BOOST_CHECK_THROW(System::FromCanonicalEncoding(""), std::runtime_error);
	BOOST_CHECK_THROW(System::FromCanonicalEncoding("not an encoding at all"), std::runtime_error);

	BOOST_CHECK_THROW(DecodeCanonicalTree("(frob 1:x)"), std::runtime_error);      // unknown node kind
	BOOST_CHECK_THROW(DecodeCanonicalTree("(sum +(var 1:x) +#5)"), std::runtime_error);   // dangling back-reference
	BOOST_CHECK_THROW(DecodeCanonicalTree("(var 7:x)"), std::runtime_error);       // netstring past the end
	BOOST_CHECK_THROW(DecodeCanonicalTree("(var 1:x) (var 1:y)"), std::runtime_error);   // two trees
}

BOOST_AUTO_TEST_SUITE_END()
