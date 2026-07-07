//This file is part of Bertini 2.
//
//bertini2/io/parsing/unicode_ident.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/unicode_ident.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/unicode_ident.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/unicode_ident.hpp

 \brief UTF-8-aware identifier building blocks for the classic (Bertini 1) Qi grammar.

 The classic parser keeps byte iterators (`std::string::const_iterator`)
 everywhere -- the symbol tables (`qi::symbols<char,...>`) and error handling are
 all keyed on `char`.  To accept Unicode *letters* (Ω, α, CJK, ...) as variable
 names without churning that whole pipeline, we make only the identifier rule
 UTF-8-aware: code points are transiently decoded to UTF-32 via
 `boost::u8_to_u32_iterator` and classified by Unicode character properties
 (never the platform locale / `wchar_t`, which is the cp1252/latin-1 footgun).
 Because the source string is already UTF-8, a matched identifier's raw bytes
 *are* the desired `std::string` attribute.

 Two self-contained Spirit primitives are provided so composition never needs a
 unary/binary Spirit operator applied to a bare custom primitive (which requires
 the terminal machinery): utf8_identifier_parser matches a whole identifier and
 exposes it as a `std::string`; utf8_ident_boundary_parser is a zero-width
 negative lookahead used to stop a known symbol from matching a prefix of a
 longer identifier.
 */

#pragma once

#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/regex/pending/unicode_iterator.hpp>

#include "bertini2/naming.hpp"   // IsIdentStart / IsIdentCont -- the shared name policy

namespace bertini {
namespace parsing {
namespace classic {

// IsIdentStart / IsIdentCont live in namespace bertini (bertini2/naming.hpp) so the
// grammar and the object model share one definition of a legal name; they resolve
// here by enclosing-namespace lookup.

/// \cond UNICODE_IDENT_PARSERS

// Match a whole UTF-8 identifier (IsIdentStart, then zero or more IsIdentCont)
// and expose the matched raw bytes as a std::string.  Self-contained: decodes
// contiguous code points itself, so it never sub-skips inside an identifier even
// when used in a rule that carries a space skipper.
struct utf8_identifier_parser
    : boost::spirit::qi::primitive_parser<utf8_identifier_parser>
{
	template <typename Context, typename Iterator>
	struct attribute { typedef std::string type; };

	template <typename Iterator, typename Context, typename Skipper, typename Attribute>
	bool parse(Iterator& first, Iterator const& last,
	           Context&, Skipper const& skipper, Attribute& attr) const
	{
		boost::spirit::qi::skip_over(first, last, skipper);
		Iterator const begin = first;
		Iterator it = first;
		try
		{
			if (it == last)
				return false;
			{
				boost::u8_to_u32_iterator<Iterator> u(it, it, last);
				if (!IsIdentStart(*u))
					return false;
				++u;
				it = u.base();
			}
			while (it != last)
			{
				boost::u8_to_u32_iterator<Iterator> u(it, it, last);
				if (!IsIdentCont(*u))
					break;
				++u;
				it = u.base();
			}
		}
		catch (...)
		{
			return false; // malformed UTF-8 -> clean parse failure
		}
		std::string matched(begin, it);
		boost::spirit::traits::assign_to(matched, attr);
		first = it;
		return true;
	}

	template <typename Context>
	boost::spirit::info what(Context&) const
	{
		return boost::spirit::info("utf8_identifier");
	}
};

// Zero-width negative lookahead: succeeds (consuming nothing beyond the skipper)
// iff the next code point is NOT an identifier-continuation code point (or the
// input has ended).  Replaces the ASCII `!qi::alnum` guard so a following
// Unicode letter (e.g. `α` after a known `Ω`) also blocks a greedy symbol match.
struct utf8_ident_boundary_parser
    : boost::spirit::qi::primitive_parser<utf8_ident_boundary_parser>
{
	template <typename Context, typename Iterator>
	struct attribute { typedef boost::spirit::unused_type type; };

	template <typename Iterator, typename Context, typename Skipper, typename Attribute>
	bool parse(Iterator& first, Iterator const& last,
	           Context&, Skipper const& skipper, Attribute&) const
	{
		boost::spirit::qi::skip_over(first, last, skipper);
		if (first == last)
			return true; // nothing follows -> boundary holds
		try
		{
			boost::u8_to_u32_iterator<Iterator> u(first, first, last);
			return !IsIdentCont(*u); // followed by an ident char -> boundary fails
		}
		catch (...)
		{
			return true; // malformed following bytes -> treat as a boundary
		}
	}

	template <typename Context>
	boost::spirit::info what(Context&) const
	{
		return boost::spirit::info("utf8_ident_boundary");
	}
};

/// \endcond

} // namespace classic
} // namespace parsing
} // namespace bertini
