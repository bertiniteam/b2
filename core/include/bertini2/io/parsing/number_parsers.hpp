//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/number_parsers.hpp
 
 \brief Provides functions for parsing numbers in bertini.
 */

#pragma once

#include "bertini2/io/parsing/qi_files.hpp"
#include "bertini2/io/parsing/number_rules.hpp"



namespace bertini{
	namespace parsers{

		namespace classic{
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, double& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				double rN = 0.0;
				bool r = phrase_parse(first, last,
									  (
									   double_[ref(rN) = _1]
									   ),
									  space);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				c = rN;
				return r;
			}
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, std::complex<double>& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				double rN = 0.0;
				double iN = 0.0;
				bool r = phrase_parse(first, last,
									  (
									   double_[ref(rN) = _1]
									   >> -(double_[ref(iN) = _1])
									   |   double_[ref(rN) = _1]
									   ),
									  space);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				c = std::complex<double>(rN, iN);
				return r;
			}
			
			
			
			
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, mpfr_float& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				MpfrFloatParser<Iterator> S;
				
				mpfr_float rN {0};
				bool r = phrase_parse(first, last,
									  S,
									  space,
									  rN
									  );
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				c = rN;
				return r;
			}
			
			
			
			
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, mpfr& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				MpfrComplexParser<Iterator> S;
				
				mpfr rN {};
				bool r = phrase_parse(first, last,
									  S,
									  space,
									  rN);
				
				std::cout << rN.real().str() << std::endl;
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				c = rN;
				return r;
			}
			
		} // re: namespace classic



		namespace cplusplus{
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, double& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				double rN = 0.0;
				bool r = phrase_parse(first, last,
									  (
									   double_[ref(rN) = _1]
									   ),
									  space);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				c = rN;
				return r;
			}
			
			
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, std::complex<double>& c)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				double rN = 0.0;
				double iN = 0.0;
				bool r = phrase_parse(first, last,
									  (
									   '(' >> double_[ref(rN) = _1]
									   >> -(',' >> double_[ref(iN) = _1]) >> ')'
									   |   double_[ref(rN) = _1]
									   ),
									  space);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				c = std::complex<double>(rN, iN);
				return r;
			}
		} // re: namespace cplusplus
	} //re: namespace parsing
} //re: namespace bertini
