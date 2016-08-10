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
\file bertini2/io/parsers.hpp 

\brief Provides the parsers for bertini2.
*/

#pragma once




#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/phoenix/object/construct.hpp>
#include <boost/bind.hpp>

#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include "bertini2/mpfr_complex.hpp"




BOOST_FUSION_ADAPT_ADT(
		    bertini::complex,
		    (obj.real(), obj.real(val))
		    (obj.imag(), obj.imag(val)))


namespace bertini{
	namespace parsers{

		

		namespace qi = ::boost::spirit::qi;
		namespace ascii = ::boost::spirit::ascii;

		/** 
		/brief Struct that holds the rules for parsing mpfr numbers
		
		Use: Include this struct with the rest of your rules.
		*/
		template<typename Iterator>
		struct MPParserRules
		{
			MPParserRules()
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
				number_with_no_point_ [_val += _1]
				>> -exponent_notation_ [_val += _1];
				
				
				long_number_string_.name("long_number_string_");
				long_number_string_ = eps[_val = std::string()] >>
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
				+(qi::char_(L'0',L'9')[_val += _1])
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
				
			}
			
			
			
			// these rules all produce strings which are fed into numbers.
			qi::rule<Iterator, std::string()> number_string_, integer_string_, long_number_string_, number_with_digits_before_point_,
			number_with_digits_after_point_, number_with_no_point_, exponent_notation_;
			
		}; //re: MPParserRules


		template<typename Iterator, typename Skipper = ascii::space_type>
		struct MpfrFloatParser : qi::grammar<Iterator, mpfr_float(), boost::spirit::ascii::space_type>
		{
			MpfrFloatParser() : MpfrFloatParser::base_type(root_rule_,"MpfrFloatParser")
			{
				namespace phx = boost::phoenix;
				using qi::_1;
				using qi::_val;

				root_rule_.name("mpfr_float");
				root_rule_ = mpfr_rules_.number_string_
								[ phx::bind( 
										[]
										(mpfr_float & B, std::string P)
										{
											auto prev_prec = DefaultPrecision();

											DefaultPrecision(P.size());
											B = mpfr_float{P};
											DefaultPrecision(prev_prec);
											assert(B.precision() == P.size());
										},
										_val,_1
									)
								];
			}

			qi::rule<Iterator, mpfr_float(), Skipper > root_rule_;
			MPParserRules<Iterator> mpfr_rules_;
		};





		struct Classic{

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



		    template<typename Iterator, typename Skipper = ascii::space_type>
			struct MpfrComplexParser : qi::grammar<Iterator, mpfr(), boost::spirit::ascii::space_type>
			{
				MpfrComplexParser() : MpfrComplexParser::base_type(root_rule_,"MpfrComplexParser")
				{
					using qi::_1;
					using qi::_2;
					using qi::_val;

					root_rule_ = 
						(mpfr_float_ >> mpfr_float_)
						[ _val = boost::phoenix::construct<bertini::complex>(_1, _2) ];
				}

				qi::rule<Iterator, mpfr(), Skipper > root_rule_;
				MpfrFloatParser<Iterator> mpfr_float_;
			};
			    

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

		        if (!r || first != last) // fail if we did not get a full match
		            return false;

		        c = rN;
		        return r;
		    }

		};



		struct CPlusPlus{

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
		};

		
	}
}

