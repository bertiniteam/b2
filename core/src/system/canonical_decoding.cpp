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

\brief The System-level canonical-encoding reader: the inverse of content_identity.cpp's
encoder, section by section (`b2sysenc/1`).  Every section the encoder writes is read
here in the same order; a section the encoder gains must be read here in the same
commit, or the round trip (digest equality) breaks loudly.
*/

#include "bertini2/system/system.hpp"

#include "bertini2/function_tree/canonical.hpp"
#include "bertini2/function_tree/canonical_decoding.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bertini {

namespace {

    using node::DecodingContext;
    using node::DecodingCursor;

    char const* OrderName(node::MonomialOrder o)
    {
        switch (o)
        {
            case node::MonomialOrder::Lex: return "Lex";
            case node::MonomialOrder::RevLex: return "RevLex";
            case node::MonomialOrder::GrevLex: return "GrevLex";
        }
        throw std::runtime_error("unencodable MonomialOrder");
    }

    VariableGroupType GroupTypeOf(std::string const& tag, DecodingCursor const& cur)
    {
        if (tag == "A") return VariableGroupType::Affine;
        if (tag == "H") return VariableGroupType::Homogeneous;
        if (tag == "U") return VariableGroupType::Ungrouped;
        cur.Fail("unknown variable group type \"" + tag + "\"");
    }

    // A group of variables, by name: `<count> <name> <name> ...` (names are netstrings).
    VariableGroup ReadVariableGroup(DecodingCursor& cur)
    {
        auto const count = cur.ReadUnsigned();
        VariableGroup group;
        group.reserve(static_cast<std::size_t>(count));
        for (unsigned long long ii = 0; ii < count; ++ii)
        {
            cur.Expect(' ');
            group.push_back(node::Variable::Make(cur.ReadNetstring()));
        }
        return group;
    }

    // One exact complex_mp entry: `<precision> <re> <im>`, rebuilt at its stored precision.
    complex_mp ReadComplexEntry(DecodingCursor& cur)
    {
        auto const precision = static_cast<unsigned>(cur.ReadUnsigned());
        cur.Expect(' ');
        auto const re = cur.ReadWord();
        cur.Expect(' ');
        auto const im = cur.ReadWord();
        real_mp const real_part(re, precision);
        real_mp const imag_part(im, precision);
        complex_mp z;
        z.precision(precision);
        z.real(real_part);
        z.imag(imag_part);
        return z;
    }

    // An exact matrix: `<rows>x<cols>` then the entries in row-major order.
    Mat<complex_mp> ReadComplexMatrix(DecodingCursor& cur)
    {
        auto const rows = cur.ReadUnsigned();
        cur.Expect('x');
        auto const cols = cur.ReadUnsigned();
        Mat<complex_mp> m(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
        for (Eigen::Index r = 0; r < m.rows(); ++r)
            for (Eigen::Index c = 0; c < m.cols(); ++c)
            {
                cur.Expect(' ');
                m(r, c) = ReadComplexEntry(cur);
            }
        return m;
    }

    // Per-row multidegree tables: `<rows> [<n> d d ...] [<n> ...] ...`
    std::vector<std::vector<int>> ReadMultidegrees(DecodingCursor& cur)
    {
        auto const rows = cur.ReadUnsigned();
        std::vector<std::vector<int>> mds;
        mds.reserve(static_cast<std::size_t>(rows));
        for (unsigned long long r = 0; r < rows; ++r)
        {
            cur.Expect(" [");
            auto const n = cur.ReadUnsigned();
            std::vector<int> row;
            row.reserve(static_cast<std::size_t>(n));
            for (unsigned long long ii = 0; ii < n; ++ii)
            {
                cur.Expect(' ');
                row.push_back(static_cast<int>(cur.ReadInt()));
            }
            cur.Expect(']');
            mds.push_back(std::move(row));
        }
        return mds;
    }

    // The path variable slot: a netstring name, or '-' for none.
    System::Var ReadOptionalVariable(DecodingCursor& cur)
    {
        if (cur.TryConsume('-'))
            return nullptr;
        return node::Variable::Make(cur.ReadNetstring());
    }

    System::NE AsNamedExpression(std::shared_ptr<node::Node> const& n, DecodingCursor const& cur)
    {
        auto named = std::dynamic_pointer_cast<node::NamedExpression>(n);
        if (!named)
            cur.Fail("expected a named expression");
        return named;
    }

} // unnamed namespace


System System::FromCanonicalEncoding(std::string const& text)
{
    DecodingCursor cursor(text);
    DecodingContext ctx;
    System sys;
    sys.DecodeCanonicalFrom(cursor, ctx);
    cursor.SkipWhitespace();
    if (!cursor.AtEnd())
        cursor.Fail("trailing bytes after the system");
    return sys;
}


void System::DecodeCanonicalFrom(DecodingCursor& cur, DecodingContext& ctx)
{
    // 1. format version + the session-global canonicalization settings the trees were shaped by
    cur.Expect("b2sysenc/");
    auto const version = cur.ReadUnsigned();
    std::string const version_token = "b2sysenc/" + std::to_string(version);
    if (version_token != SystemEncodingVersion)
        throw std::runtime_error("canonical encoding: version " + version_token
            + " is not readable by this build, which reads " + SystemEncodingVersion);
    cur.Expect(" order=");
    auto const order = cur.ReadWord();
    cur.Expect(" canon=");
    auto const canon = cur.ReadUnsigned();
    cur.Expect(" powerfold=");
    auto const powerfold = cur.ReadUnsigned();
    cur.Expect('\n');
    {
        std::string const session_order = OrderName(node::CurrentMonomialOrder());
        auto const session_canon = node::CanonicalizeByDefault() ? 1u : 0u;
        auto const session_fold = node::PowerFoldByDefault() ? 1u : 0u;
        if (order != session_order || canon != session_canon || powerfold != session_fold)
            throw std::runtime_error("canonical encoding: written under canonicalization settings order="
                + order + " canon=" + std::to_string(canon) + " powerfold=" + std::to_string(powerfold)
                + ", but this session uses order=" + session_order + " canon=" + std::to_string(session_canon)
                + " powerfold=" + std::to_string(session_fold)
                + "; decoding under other settings would not reproduce the digest");
    }

    // start from nothing: the default constructor may have seeded a polynomial block
    ungrouped_variables_.clear();
    variable_groups_.clear();
    hom_variable_groups_.clear();
    homogenizing_variables_.clear();
    pre_homogenization_functions_.clear();
    path_variable_ = nullptr;
    have_path_variable_ = false;
    implicit_parameters_.clear();
    explicit_parameters_.clear();
    blocks_.clear();
    patch_ = Patch();
    is_patched_ = false;
    time_order_of_variable_groups_.clear();
    have_ordering_ = false;
    variable_ordering_.clear();
    is_differentiated_ = false;
    sealed_digest_.reset();

    // 2. variable structure
    cur.Expect("timeorder");
    while (cur.TryConsume(' '))
        time_order_of_variable_groups_.push_back(GroupTypeOf(cur.ReadWord(), cur));
    cur.Expect('\n');

    cur.Expect("affine ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            variable_groups_.push_back(ReadVariableGroup(cur));
        }
    }
    cur.Expect('\n');

    cur.Expect("hom ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            hom_variable_groups_.push_back(ReadVariableGroup(cur));
        }
    }
    cur.Expect('\n');

    cur.Expect("ungrouped ");
    ungrouped_variables_ = ReadVariableGroup(cur);
    cur.Expect('\n');

    cur.Expect("homvars ");
    homogenizing_variables_ = ReadVariableGroup(cur);
    cur.Expect('\n');

    cur.Expect("auxgroups ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            auxiliary_variable_groups_.push_back(static_cast<unsigned>(cur.ReadUnsigned()));
        }
    }
    cur.Expect('\n');

    cur.Expect("auxcoords ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            auxiliary_coordinates_.push_back(static_cast<unsigned>(cur.ReadUnsigned()));
        }
    }
    cur.Expect('\n');

    cur.Expect("pathvar ");
    path_variable_ = ReadOptionalVariable(cur);
    have_path_variable_ = static_cast<bool>(path_variable_);
    cur.Expect('\n');

    // 3. parameters
    cur.Expect("implicit ");
    implicit_parameters_ = ReadVariableGroup(cur);
    cur.Expect('\n');

    cur.Expect("explicit ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            explicit_parameters_.push_back(AsNamedExpression(node::DecodeCanonical(cur, ctx), cur));
        }
    }
    cur.Expect('\n');

    // 4. blocks, in stored order
    cur.Expect("blocks ");
    auto const num_blocks = cur.ReadUnsigned();
    cur.Expect('\n');
    for (unsigned long long b = 0; b < num_blocks; ++b)
    {
        cur.Expect("(block ");
        auto const kind = cur.ReadWord();
        if (kind == "poly")
        {
            cur.Expect(' ');
            auto const n = cur.ReadUnsigned();
            blocks::PolynomialBlock block;
            for (unsigned long long ii = 0; ii < n; ++ii)
            {
                cur.Expect(' ');
                block.AddFunction(node::DecodeCanonical(cur, ctx));
            }
            cur.Expect(")\n");
            blocks_.emplace_back(std::move(block));
        }
        else if (kind == "linforms")
        {
            cur.Expect(' ');
            auto const num_vars = cur.ReadUnsigned();
            cur.Expect(' ');
            auto const homogenized = cur.ReadUnsigned();
            cur.Expect(' ');
            auto coefficients = ReadComplexMatrix(cur);
            cur.Expect(")\n");
            blocks_.emplace_back(blocks::LinearFormsBlock(static_cast<size_t>(num_vars),
                                                          std::move(coefficients), homogenized == 1));
        }
        else if (kind == "prodlin")
        {
            cur.Expect(' ');
            auto const num_vars = cur.ReadUnsigned();
            cur.Expect(' ');
            auto const num_factors = cur.ReadUnsigned();
            std::vector<Mat<complex_mp>> factors;
            factors.reserve(static_cast<std::size_t>(num_factors));
            for (unsigned long long ii = 0; ii < num_factors; ++ii)
            {
                cur.Expect(' ');
                factors.push_back(ReadComplexMatrix(cur));
            }
            cur.Expect(")\n");
            blocks_.emplace_back(blocks::ProductsOfLinearsBlock(static_cast<size_t>(num_vars), std::move(factors)));
        }
        else if (kind == "randomization")
        {
            cur.Expect(" homogenized=");
            auto const homogenized = cur.ReadUnsigned();
            cur.Expect(" groups=");
            auto const num_groups = cur.ReadUnsigned();
            cur.Expect(" R=");
            auto R = ReadComplexMatrix(cur);
            cur.Expect(" target=");
            auto target_mds = ReadMultidegrees(cur);
            cur.Expect(" operandmd=");
            auto operand_mds = ReadMultidegrees(cur);
            cur.Expect(" homvars=");
            std::vector<Var> hom_vars;
            {
                auto const n = cur.ReadUnsigned();
                for (unsigned long long ii = 0; ii < n; ++ii)
                {
                    cur.Expect(' ');
                    hom_vars.push_back(node::Variable::Make(cur.ReadNetstring()));
                }
            }
            cur.Expect(" operand=(\n");
            auto operand = std::make_shared<System>();
            operand->DecodeCanonicalFrom(cur, ctx);   // recursion, same context
            cur.Expect("))\n");
            blocks::RandomizationBlock<System> block(std::move(operand), std::move(R),
                                                     std::move(target_mds), std::move(operand_mds),
                                                     static_cast<size_t>(num_groups));
            if (homogenized == 1)
                block.RestoreHomogenization(std::move(hom_vars));
            blocks_.emplace_back(std::move(block));
        }
        else if (kind == "blend")
        {
            cur.Expect(" pathvar ");
            auto path_variable = ReadOptionalVariable(cur);
            cur.Expect(" coefficients ");
            auto const num_coefficients = cur.ReadUnsigned();
            std::vector<Nd> coefficients;
            coefficients.reserve(static_cast<std::size_t>(num_coefficients));
            for (unsigned long long ii = 0; ii < num_coefficients; ++ii)
            {
                cur.Expect(' ');
                coefficients.push_back(node::DecodeCanonical(cur, ctx));
            }
            cur.Expect(" operands ");
            auto const num_operands = cur.ReadUnsigned();
            std::vector<std::shared_ptr<const System>> operands;
            operands.reserve(static_cast<std::size_t>(num_operands));
            for (unsigned long long ii = 0; ii < num_operands; ++ii)
            {
                cur.Expect(" (\n");
                auto operand = std::make_shared<System>();
                operand->DecodeCanonicalFrom(cur, ctx);   // recursion, same context
                cur.Expect(')');
                operands.push_back(std::move(operand));
            }
            cur.Expect(")\n");
            blocks_.emplace_back(blocks::BlendBlock<System>(std::move(path_variable), std::move(coefficients),
                                                            std::move(operands)));
        }
        else
            cur.Fail("unknown block kind \"" + kind + "\"");
    }

    // 5. the patch
    cur.Expect("patched ");
    {
        auto const patched = cur.ReadUnsigned();
        if (patched == 1)
        {
            cur.Expect(" sizes ");
            auto const num_sizes = cur.ReadUnsigned();
            std::vector<unsigned> sizes;
            sizes.reserve(static_cast<std::size_t>(num_sizes));
            for (unsigned long long ii = 0; ii < num_sizes; ++ii)
            {
                cur.Expect(' ');
                sizes.push_back(static_cast<unsigned>(cur.ReadUnsigned()));
            }
            cur.Expect(" coefficients ");
            auto const num_vectors = cur.ReadUnsigned();
            std::vector<Vec<complex_mp>> coefficients;
            coefficients.reserve(static_cast<std::size_t>(num_vectors));
            for (unsigned long long ii = 0; ii < num_vectors; ++ii)
            {
                cur.Expect(" [");
                auto const n = cur.ReadUnsigned();
                Vec<complex_mp> vec(static_cast<Eigen::Index>(n));
                for (Eigen::Index jj = 0; jj < vec.size(); ++jj)
                {
                    cur.Expect(' ');
                    vec(jj) = ReadComplexEntry(cur);
                }
                cur.Expect(']');
                coefficients.push_back(std::move(vec));
            }
            patch_ = Patch(sizes, coefficients);
            is_patched_ = true;
        }
        else if (patched != 0)
            cur.Fail("the patched flag must be 0 or 1");
    }
    cur.Expect('\n');

    // 6. the pre-homogenization snapshot
    cur.Expect("prehom ");
    {
        auto const n = cur.ReadUnsigned();
        for (unsigned long long ii = 0; ii < n; ++ii)
        {
            cur.Expect(' ');
            pre_homogenization_functions_.push_back(node::DecodeCanonical(cur, ctx));
        }
    }
    cur.Expect('\n');
}

} // namespace bertini
