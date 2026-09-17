//This file is part of Bertini 2.
//
//canonical_decoding.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_decoding.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_decoding.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical_decoding.cpp

\brief The canonical-encoding reader for expression trees (ADR-0042): the inverse of
the encoder in canonical_encoding.cpp, node kind by node kind.  If the encoder learns a
new kind, this reader must learn it in the same commit.
*/

#include "bertini2/function_tree/canonical_decoding.hpp"

#include "bertini2/function_tree.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace bertini {
namespace node {

// ---- the cursor ----

char DecodingCursor::Peek() const
{
    if (AtEnd())
        Fail("unexpected end of the encoding");
    return text_[pos_];
}

void DecodingCursor::Expect(char c)
{
    if (AtEnd() || text_[pos_] != c)
        Fail(std::string("expected '") + c + "'");
    ++pos_;
}

void DecodingCursor::Expect(std::string_view lit)
{
    if (text_.compare(pos_, lit.size(), lit) != 0)
        Fail("expected \"" + std::string(lit) + "\"");
    pos_ += lit.size();
}

bool DecodingCursor::TryConsume(char c)
{
    if (AtEnd() || text_[pos_] != c)
        return false;
    ++pos_;
    return true;
}

void DecodingCursor::SkipWhitespace()
{
    while (!AtEnd() && (text_[pos_] == ' ' || text_[pos_] == '\n'))
        ++pos_;
}

std::string DecodingCursor::ReadWord()
{
    auto const start = pos_;
    while (!AtEnd())
    {
        char const c = text_[pos_];
        if (c == ' ' || c == '\n' || c == '(' || c == ')' || c == ']')
            break;
        ++pos_;
    }
    if (pos_ == start)
        Fail("expected a word");
    return text_.substr(start, pos_ - start);
}

std::string DecodingCursor::ReadNetstring()
{
    auto const length = ReadUnsigned();
    Expect(':');
    if (length > text_.size() - pos_)
        Fail("netstring runs past the end of the encoding");
    auto const start = pos_;
    pos_ += static_cast<std::size_t>(length);
    return text_.substr(start, static_cast<std::size_t>(length));
}

unsigned long long DecodingCursor::ReadUnsigned()
{
    auto const start = pos_;
    while (!AtEnd() && std::isdigit(static_cast<unsigned char>(text_[pos_])))
        ++pos_;
    if (pos_ == start)
        Fail("expected a non-negative integer");
    return std::stoull(text_.substr(start, pos_ - start));
}

long long DecodingCursor::ReadInt()
{
    auto const start = pos_;
    if (!AtEnd() && (text_[pos_] == '-' || text_[pos_] == '+'))
        ++pos_;
    while (!AtEnd() && std::isdigit(static_cast<unsigned char>(text_[pos_])))
        ++pos_;
    if (pos_ == start || (pos_ == start + 1 && !std::isdigit(static_cast<unsigned char>(text_[start]))))
        Fail("expected an integer");
    return std::stoll(text_.substr(start, pos_ - start));
}

void DecodingCursor::Fail(std::string const& what) const
{
    auto const from = pos_ > 20 ? pos_ - 20 : 0;
    auto const snippet = text_.substr(from, 60);
    std::ostringstream msg;
    msg << "canonical encoding: " << what << " at byte " << pos_
        << " (near \"" << snippet << "\")";
    throw std::runtime_error(msg.str());
}


// ---- the node reader ----

namespace {

    using Nd = std::shared_ptr<Node>;

    // A complex_mp literal at its STORED precision: the precision is set first so the
    // digits are taken exactly as written, never rounded through the session default.
    complex_mp ComplexAtPrecision(unsigned precision, std::string const& re, std::string const& im)
    {
        real_mp const real_part(re, precision);
        real_mp const imag_part(im, precision);
        complex_mp z;
        z.precision(precision);
        z.real(real_part);
        z.imag(imag_part);
        return z;
    }

    class NodeReader
    {
    public:
        NodeReader(DecodingCursor& cursor, DecodingContext& ctx) : cur_(cursor), ctx_(ctx) {}

        Nd Read()
        {
            cur_.SkipWhitespace();
            if (cur_.TryConsume('#'))
            {
                auto const index = cur_.ReadUnsigned();
                if (index >= ctx_.by_index.size() || !ctx_.by_index[static_cast<std::size_t>(index)])
                    cur_.Fail("back-reference #" + std::to_string(index) + " names a node not yet read");
                return ctx_.by_index[static_cast<std::size_t>(index)];
            }

            cur_.Expect('(');
            auto const tag = cur_.ReadWord();

            // the node's index is claimed before its children are read, as the encoder
            // numbered it before descending
            auto const my_index = ctx_.by_index.size();
            ctx_.by_index.push_back(nullptr);

            Nd result = ReadBody(tag);
            ctx_.by_index[my_index] = result;
            return result;
        }

    private:
        Nd ReadBody(std::string const& tag)
        {
            if (tag == "var")
            {
                cur_.Expect(' ');
                auto const name = cur_.ReadNetstring();
                cur_.Expect(')');
                return Variable::Make(name);
            }
            if (tag == "diff")
            {
                cur_.Expect(' ');
                auto const name = cur_.ReadNetstring();
                cur_.Expect(')');
                // a differential is named after its variable (Variable::Differentiate does the same)
                return Differential::Make(Variable::Make(name), name);
            }
            if (tag == "int")
            {
                cur_.Expect(' ');
                auto const digits = cur_.ReadWord();
                cur_.Expect(')');
                return Integer::Make(mpz_int(digits));
            }
            if (tag == "rat")
            {
                cur_.Expect(' ');
                auto const re = cur_.ReadWord();
                cur_.Expect(' ');
                auto const im = cur_.ReadWord();
                cur_.Expect(')');
                return Rational::Make(mpq_rational(re), mpq_rational(im));
            }
            if (tag == "cplx")
            {
                cur_.Expect(' ');
                auto const precision = cur_.ReadUnsigned();
                cur_.Expect(' ');
                auto const re = cur_.ReadWord();
                cur_.Expect(' ');
                auto const im = cur_.ReadWord();
                cur_.Expect(')');
                return Complex::Make(ComplexAtPrecision(static_cast<unsigned>(precision), re, im));
            }
            if (tag == "pi")
            {
                cur_.Expect(')');
                return Pi();
            }
            if (tag == "e")
            {
                cur_.Expect(')');
                return E();
            }
            if (tag == "named")
            {
                cur_.Expect(' ');
                auto const name = cur_.ReadNetstring();
                cur_.Expect(' ');
                auto entry = Read();
                cur_.Expect(')');
                return NamedExpression::Make(entry, name);
            }
            if (tag == "sum")
                return ReadNary('+', '-', /*is_sum=*/true);
            if (tag == "mul")
                return ReadNary('*', '/', /*is_sum=*/false);
            if (tag == "pow")
            {
                cur_.Expect(' ');
                auto base = Read();
                cur_.Expect(' ');
                auto exponent = Read();
                cur_.Expect(')');
                return PowerOperator::Make(base, exponent);
            }
            if (tag == "ipow")
            {
                cur_.Expect(' ');
                auto const exponent = cur_.ReadInt();
                cur_.Expect(' ');
                auto operand = Read();
                cur_.Expect(')');
                return IntegerPowerOperator::Make(operand, static_cast<int>(exponent));
            }

            // the unary family
            cur_.Expect(' ');
            auto operand = Read();
            cur_.Expect(')');
            if (tag == "neg")  return NegateOperator::Make(operand);
            if (tag == "sqrt") return SqrtOperator::Make(operand);
            if (tag == "exp")  return ExpOperator::Make(operand);
            if (tag == "log")  return LogOperator::Make(operand);
            if (tag == "sin")  return SinOperator::Make(operand);
            if (tag == "asin") return ArcSinOperator::Make(operand);
            if (tag == "cos")  return CosOperator::Make(operand);
            if (tag == "acos") return ArcCosOperator::Make(operand);
            if (tag == "tan")  return TanOperator::Make(operand);
            if (tag == "atan") return ArcTanOperator::Make(operand);

            cur_.Fail("unknown node kind \"" + tag + "\"");
        }

        // (sum +a -b ...) / (mul *a /b ...): a sign or mult-or-div mark before each operand
        Nd ReadNary(char positive, char negative, bool is_sum)
        {
            std::vector<std::pair<Nd, bool>> operands;
            while (cur_.TryConsume(' '))
            {
                char const mark = cur_.Peek();
                if (mark != positive && mark != negative)
                    cur_.Fail(std::string("expected '") + positive + "' or '" + negative + "' before an operand");
                cur_.Expect(mark);
                auto operand = Read();
                operands.emplace_back(std::move(operand), mark == positive);
            }
            cur_.Expect(')');
            if (is_sum)
                return SumOperator::Make(operands);
            return MultOperator::Make(operands);
        }

        DecodingCursor& cur_;
        DecodingContext& ctx_;
    };

} // unnamed namespace


std::shared_ptr<Node> DecodeCanonical(DecodingCursor& cursor, DecodingContext& ctx)
{
    NodeReader reader(cursor, ctx);
    return reader.Read();
}

std::shared_ptr<Node> DecodeCanonicalTree(std::string const& text)
{
    DecodingCursor cursor(text);
    DecodingContext ctx;
    auto root = DecodeCanonical(cursor, ctx);
    cursor.SkipWhitespace();
    if (!cursor.AtEnd())
        cursor.Fail("trailing bytes after the tree");
    return root;
}

} // namespace node
} // namespace bertini
