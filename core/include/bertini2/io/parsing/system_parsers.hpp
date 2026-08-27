//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/system_parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/system_parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/system_parsers.hpp
 
 \brief Provides the parsers for systems in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/system_rules.hpp"



namespace bertini {
	namespace parsing {
		
		namespace classic {
		
		
		
			/// \brief Parse a polynomial system from the iterator range into a System.
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, System& sys)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				SystemParser<Iterator> S;

				System s{};
				bool r = phrase_parse(first, last,
									  S,
									  space,
									  s);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;

				// every definition has now filled its function box: emit the bare entries.
				S.EmitDeclaredFunctions(s);
				sys = s;
				return r;
			}
		
		} // re: namespace classic
		
	}// re: namespace parsing

	inline
	System::System(std::string const& input)
	{
		System sys;

		parsing::classic::SystemParser<std::string::const_iterator> S;

		// Treat the input as UTF-8; drop a leading BOM so it is not seen as a
		// stray leading character by the grammar.
		std::string cleaned = input;
		parsing::classic::StripUTF8BOM(cleaned);
		// accept a full Bertini 1 classic file (CONFIG/INPUT wrappers), not just the bare
		// INPUT-section body that the grammar reads -- see #396
		parsing::classic::StripClassicFileWrappers(cleaned);

		std::string::const_iterator iter = cleaned.begin();
		std::string::const_iterator end = cleaned.end();
		
		bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);
		
		if (!s || iter!=end)
		{
			std::string remaining(iter, end);
			if (remaining.size() > 60)
				remaining = remaining.substr(0, 60) + "...";
			if (remaining.empty())
				remaining = "<end of input>";
			throw std::runtime_error(
				"[SystemParser] parser did not consume entire input; "
				"unparsed remainder: \"" + remaining + "\"");
		}
		
		// every definition has now filled its function box: emit the bare entries.
		S.EmitDeclaredFunctions(sys);

		using std::swap;
		swap(sys,*this);
	}
	
}// re: namespace bertini


