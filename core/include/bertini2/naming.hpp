//This file is part of Bertini 2.
//
//bertini2/naming.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/naming.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/naming.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/naming.hpp

 \brief The single policy for what makes a legal identifier (variable/symbol name).

 A name is UTF-8: it starts with a Unicode letter (ASCII `A-Z a-z`, or Ω, α, CJK,
 ...) and continues with Unicode alphanumerics or `[ ] _`.  Classification decodes
 code points to UTF-32 (`boost::u8_to_u32_iterator`) and tests Unicode character
 properties -- never `wchar_t` / the platform locale (the cp1252/latin-1 footgun).

 These predicates are shared by the classic Qi grammar (see
 io/parsing/unicode_ident.hpp) and the object model (`node::Variable`), so the
 parser and direct construction agree on exactly what a name may be.  Emoji are
 symbol-category code points, not letters, so they are rejected by default; the
 two predicates below are the single opt-in point if that ever changes.
 */

#pragma once

#include <string>
#include <stdexcept>

#include <boost/regex/pending/unicode_iterator.hpp>
#include <boost/spirit/home/support/char_encoding/unicode.hpp>

namespace bertini {

/// \brief True if code point \p cp may START an identifier: any Unicode letter
///        (ASCII `A-Z a-z`, plus Ω, α, CJK, ...).
inline bool IsIdentStart(char32_t cp)
{
	return boost::spirit::char_encoding::unicode::isalpha(cp);
}

/// \brief True if code point \p cp may CONTINUE an identifier: any Unicode
///        alphanumeric, or one of the legacy continuation characters `[ ] _`.
inline bool IsIdentCont(char32_t cp)
{
	return boost::spirit::char_encoding::unicode::isalnum(cp)
	    || cp == U'[' || cp == U']' || cp == U'_';
}

/// \brief True if \p s is a well-formed identifier: nonempty, a valid UTF-8
///        IsIdentStart code point followed by zero or more IsIdentCont code
///        points.  Rejects the empty string, leading digits, operators,
///        whitespace, and malformed UTF-8.
inline bool IsValidVariableName(std::string const& s)
{
	if (s.empty())
		return false;
	try
	{
		boost::u8_to_u32_iterator<std::string::const_iterator> it(s.begin(), s.begin(), s.end());
		boost::u8_to_u32_iterator<std::string::const_iterator> end(s.end(), s.begin(), s.end());
		if (it == end || !IsIdentStart(*it))
			return false;
		for (++it; it != end; ++it)
			if (!IsIdentCont(*it))
				return false;
		return true;
	}
	catch (...)
	{
		return false; // malformed UTF-8
	}
}

/// \brief Throw std::runtime_error if \p s is not a valid variable name (\see
///        IsValidVariableName).  A no-op for a valid name.
inline void ThrowIfInvalidVariableName(std::string const& s)
{
	if (!IsValidVariableName(s))
		throw std::runtime_error(
			"invalid variable name \"" + s + "\": a name must be a nonempty identifier -- "
			"it starts with a letter (ASCII or Unicode, e.g. x or \xCE\xA9) and continues with "
			"letters, digits, or [ ] _ .  Expressions, operators, whitespace, and leading "
			"digits are not allowed.");
}

} // namespace bertini
