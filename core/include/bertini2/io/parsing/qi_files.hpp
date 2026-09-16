//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/qi_files.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/qi_files.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/qi_files.hpp
 
 \brief Provides all the include files needed to develop a parsing file that uses boost qi.
 */

#pragma once

/// \brief Select Boost.Phoenix V3 for the Spirit grammars in this header.
#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/phoenix/object/construct.hpp>
#include <boost/bind/bind.hpp>

#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cstring>
#include <string>

#include "bertini2/io/parsing/unicode_ident.hpp"

namespace bertini {
namespace parsing {
namespace classic {

/// \brief Strip a leading UTF-8 byte-order mark (`EF BB BF`) from \p s in place,
///        if present.  Editors (notably on Windows) may prepend a BOM; the Qi
///        grammar treats input as raw UTF-8, so the BOM must be removed at the
///        read boundary or it would appear as a stray leading "character".
inline void StripUTF8BOM(std::string& s)
{
	if (s.size() >= 3 &&
	    static_cast<unsigned char>(s[0]) == 0xEF &&
	    static_cast<unsigned char>(s[1]) == 0xBB &&
	    static_cast<unsigned char>(s[2]) == 0xBF)
	{
		s.erase(0, 3);
	}
}

/// \brief Unwrap a Bertini 1 classic input FILE down to the declarations the grammar reads.
///
/// `System::to_classic_input()` emits a complete Bertini 1 file --
/// `CONFIG ... END;` then `INPUT ... END;` -- because its purpose is running the same
/// problem in Bertini 1.  The grammar reads only the INPUT-section BODY.  The two therefore
/// were not inverses, and the mismatch did not announce itself: the wrapped text matched
/// nothing and the caller received an empty System.  See #396.
///
/// Comply rather than refuse: handing a system parser a classic input file means one thing.
/// The CONFIG section holds tracking settings a System does not carry, so it is dropped.
/// Text that is already a bare body is returned untouched, so this is a no-op for every
/// caller that was working before.
inline void StripClassicFileWrappers(std::string& s)
{
	auto skip_ws = [](std::string const& t, size_t i) {
		while (i < t.size() && std::isspace(static_cast<unsigned char>(t[i]))) ++i;
		return i;
	};
	// case-insensitive check for `word` at position i, followed by a non-identifier char
	auto keyword_at = [](std::string const& t, size_t i, char const* word) {
		size_t n = std::strlen(word);
		if (i + n > t.size()) return false;
		for (size_t k = 0; k < n; ++k)
			if (std::toupper(static_cast<unsigned char>(t[i+k])) != word[k]) return false;
		if (i + n < t.size())
		{
			unsigned char c = static_cast<unsigned char>(t[i+n]);
			if (std::isalnum(c) || c=='_') return false;
		}
		return true;
	};

	size_t i = skip_ws(s, 0);

	// a leading CONFIG section, if present, runs to its first END;
	if (keyword_at(s, i, "CONFIG"))
	{
		size_t e = s.find("END;", i);
		if (e == std::string::npos)
			return;                       // malformed; let the grammar report it
		i = skip_ws(s, e + 4);
	}

	// an INPUT section wrapper, if present: drop the keyword and the matching trailing END;
	if (keyword_at(s, i, "INPUT"))
	{
		size_t body = skip_ws(s, i + 5);
		size_t e = s.rfind("END;");
		if (e == std::string::npos || e < body)
			return;                       // malformed; let the grammar report it
		s = s.substr(body, e - body);
		return;
	}

	if (i != 0)
		s = s.substr(i);
}

/// \brief Format a human-readable parse-error message (line, column, expected, and found text).
inline std::string FormatParseError(
    std::string::const_iterator begin,
    std::string::const_iterator end,
    std::string::const_iterator err_pos,
    boost::spirit::info const& what,
    std::string const& parser_name)
{
    int line = 1 + (int)std::count(begin, err_pos, '\n');
    auto rev = std::find(std::make_reverse_iterator(err_pos),
                         std::make_reverse_iterator(begin), '\n');
    int col = 1 + (int)std::distance(rev.base(), err_pos);

    std::string remaining(err_pos, end);
    if (remaining.size() > 60)
        remaining = remaining.substr(0, 60) + "...";
    if (remaining.empty())
        remaining = "<end of input>";

    std::ostringstream oss;
    oss << "[" << parser_name << "] parse error at line " << line
        << ", col " << col << ":\n"
        << "  expected: " << what << "\n"
        << "  found:    \"" << remaining << "\"\n";
    return oss.str();
}

/// \brief Report a parse error (via FormatParseError) to the log, for use in a Spirit on_error handler.
inline void ReportParseError(
    std::string::const_iterator begin,
    std::string::const_iterator end,
    std::string::const_iterator err_pos,
    boost::spirit::info const& what,
    std::string const& parser_name)
{
    auto msg = FormatParseError(begin, end, err_pos, what, parser_name);
    std::cerr << msg;
    throw std::runtime_error(msg);
}

} // namespace classic
} // namespace parsing
} // namespace bertini
