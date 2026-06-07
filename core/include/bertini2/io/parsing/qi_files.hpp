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
#include <string>

namespace bertini {
namespace parsing {
namespace classic {

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
