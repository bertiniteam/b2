//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/number_rules.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/number_rules.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/io/parsing/number_rules.hpp 

\brief Provides the parsersing rules for numbers in bertini2.
*/

#pragma once




#include "bertini2/io/parsing/qi_files.hpp"

#include "bertini2/mpfr_complex.hpp"

#include "bertini2/num_traits.hpp"


BOOST_FUSION_ADAPT_ADT(
		    bertini::complex,
		    (obj.real(), obj.real(val))
		    (obj.imag(), obj.imag(val)))


namespace bertini{
	namespace parsing{

		

		namespace qi = ::boost::spirit::qi;
		namespace ascii = ::boost::spirit::ascii;

		namespace rules{
		/** 
		/brief Struct that holds the rules for parsing mpfr numbers
		
		Use: Include this struct with the rest of your rules.
		*/
		template<typename Iterator>
		struct LongNum
		{
			LongNum()
			{
				namespace phx = boost::phoenix;
				using qi::_1;
				using qi::_2;
				using qi::_3;
				using qi::_4;
				using qi::_val;
				using qi::eps;
				using qi::lit;
				using qi::char_;
				using qi::omit;
				using boost::spirit::lexeme;
				using boost::spirit::as_string;
				using boost::spirit::ascii::no_case;
				
				number_string_.name("number_");
				number_string_ =
				long_number_string_
				|
				integer_string_;
				
				integer_string_.name("integer_string_");
				integer_string_ = eps[_val = std::string()] >>
				-(qi::lit('-'))[_val += "-"] >>
				number_with_no_point_ [_val += _1]
				>> -exponent_notation_ [_val += _1];
				
				
				long_number_string_.name("long_number_string_");
				long_number_string_ = eps[_val = std::string()] >>
				-(qi::lit('-'))[_val += "-"] >>
				(
				 // 1. Read possible numbers before decimal, with possible negative
				 (number_with_digits_after_point_ [_val += _1]
				  |
				  number_with_digits_before_point_ [_val += _1] )
				 >>   // reminder -- the - before the exponent_notation here means optional
				 -exponent_notation_ [_val+=_1]// Possible scientific notation, with possible negative in exponent.
				 );
				
				
				
				number_with_digits_after_point_.name("number_with_digits_after_point_");
				number_with_digits_after_point_ = eps[_val = std::string()]
				>>
				*(qi::char_(L'0',L'9')[_val += _1])
				>>
				qi::lit('.')[_val += "."] // find a decimal point
				>>
				+(qi::char_(L'0',L'9')[_val += _1]) // find at least one digit after the point
				;
				
				
				
				
				
				number_with_digits_before_point_.name("number_with_digits_before_point_");
				number_with_digits_before_point_ = eps[_val = std::string()]
				>>
				+(qi::char_(L'0',L'9')[_val += _1])
				>>
				qi::lit('.')[_val += "."] // find a decimal point
				>>
				*(qi::char_(L'0',L'9')[_val += _1]) // find any number of digits after the point
				;
				
				
				number_with_no_point_.name("number_with_no_point_");
				number_with_no_point_ = eps[_val = std::string()]
				>>
				(-(qi::lit('-')[_val += "-"]) // then an optional minus sign
				>>
				+(qi::char_(L'0',L'9'))[_val += _1])
				;
				
				
				
				
				exponent_notation_.name("exponent_notation_");
				exponent_notation_ = eps[_val = std::string()]
				>> ( // start what the rule actually does
					(
					 qi::lit('e')[_val += "e"] // get an opening 'e'
					 |
					 qi::lit('E')[_val += "e"] // get an opening 'e'
					 )
					>>
					-(qi::lit('-')[_val += "-"]) // then an optional minus sign
					>>
					+(qi::char_(L'0',L'9')[_val += _1]) // then at least one number
					); // finish the rule off
				
				rational.name("long_rational");
				rational =
				(number_string_
				|
				(number_string_ >> char_('/') >> number_string_));

			}
			
			
			
			// these rules all produce strings which are fed into numbers.
			qi::rule<Iterator, std::string()> number_string_, integer_string_, long_number_string_, number_with_digits_before_point_,
			number_with_digits_after_point_, number_with_no_point_, exponent_notation_, rational;
			

		}; //re: LongNum


		} // namespace rules

		

		
	}
}

