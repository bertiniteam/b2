//This file is part of Bertini 2.
//
//src/function_tree/symbols/number.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/number.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/number.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/function_tree/symbols/number.hpp"

#include <sstream>
#include <string>




namespace bertini{
    namespace node{
        using ::pow;

std::shared_ptr<Node> Number::Differentiate(std::shared_ptr<Variable> const& /*v*/) const
{
    return Integer::Make(0);
}

///////////////////
//
//  INTEGERS
//
////////////////////////


void Integer::print(std::ostream & target) const
{
    target << true_value_;
}


/////////////
//
//  Floats
//
////////////////


namespace {

    // str(0) prints EVERY digit the stored value carries: streaming the number directly
    // would use the ostream's default precision (6 significant digits), silently truncating
    // every printed system -- a coefficient must round-trip exactly.
    std::string AllDigits(real_mp const& v)
    {
        return v.str(0);
    }

    // A multiprecision real in the exact Python dialect: bare when a Python literal holds it
    // exactly (an integer, or a value a double represents exactly, such as 0.5), otherwise
    // real_mp('digits', precision) so eval rebuilds it at the stored precision.
    std::string ExactPythonReal(real_mp const& v)
    {
        if (real_mp(static_cast<double>(v)) == v)
            return AllDigits(v);
        return "real_mp('" + AllDigits(v) + "', " + std::to_string(v.precision()) + ")";
    }

    // (re+im*I) or (re-im*I): self-delimiting, and the only complex spelling Bertini 1 reads.
    void PrintComplexPair(std::ostream& target, std::string const& re, std::string const& im, bool im_negative)
    {
        target << "(" << re << (im_negative ? "-" : "+") << im << "*I)";
    }

} // unnamed namespace


void Complex::print(std::ostream & target) const
{
    auto const dialect = DialectOf(target);
    auto const& re = highest_precision_value_.real();
    auto const& im = highest_precision_value_.imag();

    // real-valued floats print bare; the complex form is reserved for genuinely complex values
    if (im == 0)
    {
        target << (dialect == PrintDialect::PythonExact ? ExactPythonReal(re) : AllDigits(re));
        return;
    }

    if (dialect == PrintDialect::PythonExact)
    {
        // one node, not an expression on I: Complex('re', 'im', precision) rebuilds this very
        // leaf at its stored precision, where real_mp(re) + real_mp(im)*I would build a tree
        target << "Complex('" << AllDigits(re) << "', '" << AllDigits(im) << "', "
               << highest_precision_value_.precision() << ")";
        return;
    }
    PrintComplexPair(target, AllDigits(re), AllDigits(im < 0 ? real_mp(-im) : im), im < 0);
}


//
//  Rational
//



void Rational::print(std::ostream & target) const
{
    auto const dialect = DialectOf(target);
    auto rational_text = [&](mpq_rational const& q) -> std::string
    {
        std::ostringstream s;
        s << q;
        if (dialect == PrintDialect::PythonExact && boost::multiprecision::denominator(q) != 1)
            return "Rational('" + s.str() + "')";   // 1/3 is a double in Python; Rational('1/3') is exact
        return s.str();
    };

    // real-valued rationals print bare (p/q, textually a division -- see the precedence they
    // report); the complex form is reserved for genuinely complex values
    if (true_value_imag_ == 0)
    {
        target << rational_text(true_value_real_);
        return;
    }
    bool const im_negative = true_value_imag_ < 0;
    mpq_rational const im_abs = im_negative ? mpq_rational(-true_value_imag_) : true_value_imag_;
    PrintComplexPair(target, rational_text(true_value_real_), rational_text(im_abs), im_negative);
}


    } // re: namespace node
} // re: namespace bertini
