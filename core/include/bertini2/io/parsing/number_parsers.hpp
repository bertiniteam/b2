//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/number_parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/number_parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/number_parsers.hpp
 
 \brief Provides parsers for numbers in bertini.
 */

#pragma once

#include "bertini2/io/parsing/number_rules.hpp"



namespace bertini{
	namespace parsing{

		namespace classic{
			template<typename Iterator, typename Skipper = ascii::space_type>
			struct MpfrFloatParser : qi::grammar<Iterator, mpfr_float(), boost::spirit::ascii::space_type>
			{
				MpfrFloatParser() : MpfrFloatParser::base_type(root_rule_,"MpfrFloatParser")
				{
					using std::max;
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;

					root_rule_.name("mpfr_float");
					root_rule_ = mpfr_rules_.number_string_
									[ phx::bind( 
											[]
											(mpfr_float & B, std::string const& P)
											{
												using std::max;
												auto prev_prec = DefaultPrecision();
												auto asdf = max(prev_prec,LowestMultiplePrecision());
												auto digits = max(P.size(),static_cast<decltype(P.size())>(asdf));
												DefaultPrecision(digits);
												B = mpfr_float{P};
												DefaultPrecision(prev_prec);
												assert(B.precision() == digits);
											},
											_val,_1
										)
									];

					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("mpfr_float parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
				}

				qi::rule<Iterator, mpfr_float(), Skipper > root_rule_;
				rules::LongNum<Iterator> mpfr_rules_;
			};


			template<typename Iterator, typename Skipper = ascii::space_type>
			struct MpfrComplexParser : qi::grammar<Iterator, mpfr(), boost::spirit::ascii::space_type>
			{
				MpfrComplexParser() : MpfrComplexParser::base_type(root_rule_,"MpfrComplexParser")
				{
					using std::max;
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_val;
					
					root_rule_ =
					(mpfr_float_ >> mpfr_float_)
					[ phx::bind(
								[]
								(mpfr & B, mpfr_float const& P, mpfr_float const& Q)
								{
									auto prev_prec = DefaultPrecision();
									auto digits = max(P.precision(),Q.precision());
									
									DefaultPrecision(digits);
									B.real(P);
									B.imag(Q);
									B.precision(digits);
									DefaultPrecision(prev_prec);
									assert(B.precision() == digits);
								},
								_val,_1,_2
								)
					 ];
				}
				
				qi::rule<Iterator, mpfr(), Skipper > root_rule_;
				MpfrFloatParser<Iterator> mpfr_float_;
			};

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
