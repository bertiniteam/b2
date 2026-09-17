//This file is part of Bertini 2.
//
//print_dialect_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//print_dialect_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with print_dialect_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file print_dialect_test.cpp

\brief Printing is a family, one member per audience (ADR-0059): Bertini 1's spelling
round-trips exactly through the classic parser for every constant kind; the Python
spellings use `**` and, in the exact form, constant spellings `eval` rebuilds; a stream with
no dialect set prints Classic, as every pre-existing caller expected.
*/

#include <boost/test/unit_test.hpp>

#include <sstream>
#include <string>

#include "bertini2/function_tree.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

using namespace bertini;
using namespace bertini::node;
using Nd = std::shared_ptr<Node>;

namespace {

// A one-variable system holding `f`, and the same system rebuilt from f's classic spelling.
// The comparison is on VALUES at a multiprecision point.  Multiprecision constants come back
// bit-identical (all their digits are printed); a rational prints as p/q, which the classic
// grammar reads as a division evaluated at the working precision (Bertini 1 has no exact
// rationals), so the tolerance is the working precision's, not zero.
void CheckClassicRoundTrip(Nd const& f, std::shared_ptr<Variable> const& x)
{
    System built;
    built.AddVariableGroup(VariableGroup{x});
    built.AddFunction(f);

    std::string const text = "variable_group " + x->name() + "; function f; f = " + PrintClassic(*f) + ";";
    System parsed;
    BOOST_REQUIRE_MESSAGE(parsing::classic::parse(text.begin(), text.end(), parsed), "did not parse: " << text);

    Vec<complex_mp> point(1);
    point << complex_mp("0.37", "-1.25");
    auto const a = built.Eval(point);
    auto const b = parsed.Eval(point);
    BOOST_REQUIRE_EQUAL(a.size(), 1);
    real_mp const scale = std::max(real_mp(1), abs(a(0)));
    BOOST_CHECK_MESSAGE(abs(a(0) - b(0)) <= real_mp("1e-27") * scale,
                        "value drift through classic text: " << text << " : " << abs(a(0) - b(0)));
}

} // unnamed namespace


BOOST_AUTO_TEST_SUITE(print_dialect)

BOOST_AUTO_TEST_CASE(a_stream_with_no_dialect_prints_classic)
{
    auto x = Variable::Make("x");
    std::ostringstream out;
    out << *pow(x, 2);
    BOOST_CHECK_EQUAL(out.str(), "x^2");
    BOOST_CHECK_EQUAL(DialectOf(out) == PrintDialect::Classic, true);
}

BOOST_AUTO_TEST_CASE(powers_spell_per_dialect)
{
    auto x = Variable::Make("x"), y = Variable::Make("y");
    BOOST_CHECK_EQUAL(PrintClassic(*pow(x, 2)), "x^2");
    BOOST_CHECK_EQUAL(PrintPython(*pow(x, 2), false), "x**2");
    BOOST_CHECK_EQUAL(PrintPython(*pow(x, 2), true), "x**2");
    BOOST_CHECK_EQUAL(PrintClassic(*pow(x, -2)), "x^(-2)");
    BOOST_CHECK_EQUAL(PrintPython(*pow(x, -2), false), "x**(-2)");
    BOOST_CHECK_EQUAL(PrintClassic(*pow(x + y, y)), "(x+y)^y");
    BOOST_CHECK_EQUAL(PrintPython(*pow(x + y, y), false), "(x+y)**y");
}

BOOST_AUTO_TEST_CASE(the_dialect_rides_on_the_stream)
{
    auto x = Variable::Make("x");
    std::ostringstream out;
    SetDialect(out, PrintDialect::PythonReadable);
    out << *pow(x, 3) << " and " << *pow(x, 4);
    BOOST_CHECK_EQUAL(out.str(), "x**3 and x**4");
}

BOOST_AUTO_TEST_CASE(complex_constants_are_re_plus_im_times_I_never_a_pair)
{
    auto x = Variable::Make("x");
    auto c = Rational::Make(1, 3, 1, 2);                   // 1/3 + (1/2) i
    BOOST_CHECK_EQUAL(PrintClassic(*(c * x)), "(1/3+1/2*I)*x");
    BOOST_CHECK_EQUAL(PrintPython(*(c * x), false), "(1/3+1/2*I)*x");
    auto d = Rational::Make(1, 3, -1, 2);                  // 1/3 - (1/2) i
    BOOST_CHECK_EQUAL(PrintClassic(*(d * x)), "(1/3-1/2*I)*x");

    auto z = Complex::Make("0.5", "-0.25");
    auto const classic = PrintClassic(*z);
    BOOST_CHECK_EQUAL(classic.front(), '(');
    BOOST_CHECK(classic.find("*I)") != std::string::npos);
    BOOST_CHECK(classic.find(',') == std::string::npos);
    BOOST_CHECK(classic.find("-") != std::string::npos);   // the sign of the imaginary part
}

BOOST_AUTO_TEST_CASE(exact_python_spellings)
{
    DefaultPrecision(30);
    auto x = Variable::Make("x");
    BOOST_CHECK_EQUAL(PrintPython(*Integer::Make(-7), true), "-7");
    BOOST_CHECK_EQUAL(PrintPython(*Rational::Make(1, 3, 0, 1), true), "Rational('1/3')");
    BOOST_CHECK_EQUAL(PrintPython(*Rational::Make(1, 3, 0, 1), false), "1/3");
    BOOST_CHECK_EQUAL(PrintPython(*Rational::Make(4, 1, 0, 1), true), "4");           // integer-valued: bare
    BOOST_CHECK_EQUAL(PrintPython(*Rational::Make(1, 3, 1, 2), true), "(Rational('1/3')+Rational('1/2')*I)");
    BOOST_CHECK_EQUAL(PrintPython(*Complex::Make("0.5", "0"), true), "0.5");         // a double holds it exactly: bare
    BOOST_CHECK_EQUAL(PrintPython(*Complex::Make("2", "0"), true), "2");

    auto third = Complex::Make(complex_mp(real_mp(1) / real_mp(3)));
    auto const exact = PrintPython(*third, true);
    BOOST_CHECK_EQUAL(exact.rfind("real_mp('", 0), 0u);
    BOOST_CHECK(exact.find("', 30)") != std::string::npos);
    BOOST_CHECK_EQUAL(PrintPython(*third, false).rfind("0.3333", 0), 0u);           // readable: the digits

    auto zc = Complex::Make("0.1", "0.2");
    auto const exact_c = PrintPython(*zc, true);
    BOOST_CHECK_EQUAL(exact_c.rfind("Complex('", 0), 0u);
    BOOST_CHECK(exact_c.find(", 30)") != std::string::npos);
    BOOST_CHECK_EQUAL(PrintPython(*(zc * x), false).find("*I)*x") != std::string::npos, true);
}

BOOST_AUTO_TEST_CASE(python_spellings_never_show_a_caret_or_a_pair)
{
    auto x = Variable::Make("x"), y = Variable::Make("y");
    Nd const e = pow(x, 3) * Complex::Make("0.1", "0.2") + Rational::Make(1, 3, 1, 2) * pow(y, 2) - sqrt(x);
    for (bool exact : {false, true})
    {
        auto const s = PrintPython(*e, exact);
        BOOST_CHECK_MESSAGE(s.find('^') == std::string::npos, s);
        BOOST_CHECK_MESSAGE(s.find(',') == std::string::npos || exact, s);   // exact spellings contain ', precision)'
        BOOST_CHECK_MESSAGE(s.find("**") != std::string::npos, s);
    }
}

BOOST_AUTO_TEST_CASE(classic_spelling_round_trips_every_constant_kind_exactly)
{
    DefaultPrecision(30);
    auto x = Variable::Make("x");
    CheckClassicRoundTrip(Integer::Make(-7) * x + 3, x);
    CheckClassicRoundTrip(Rational::Make(1, 3, 0, 1) * x, x);
    CheckClassicRoundTrip(Rational::Make(1, 3, -2, 7) * pow(x, 2), x);
    CheckClassicRoundTrip(Complex::Make(complex_mp(real_mp(1) / real_mp(3))) * x, x);
    CheckClassicRoundTrip(Complex::Make(complex_mp(real_mp(1) / real_mp(3), -real_mp(2) / real_mp(7))) * x + 1, x);
    CheckClassicRoundTrip(Complex::Make("1.5e-3", "-2.0e2") * pow(x, 3), x);
    CheckClassicRoundTrip(sqrt(x) * Complex::Make("0.25", "0.75") - sin(x) / Rational::Make(1, 3, 0, 1), x);
    CheckClassicRoundTrip(pow(x + Complex::Make("0", "1"), 2), x);      // a purely imaginary constant
}

BOOST_AUTO_TEST_SUITE_END()
