//This file is part of Bertini 2.0.
//
//function_parsing.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_parsing.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_parsing.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// function_parsing.hpp:  This file contains all the grammars used to parse an input file.
//          Variable, Function
//
//      TODO:(JBC) Make these Parser classes/structs consistent.  Either they all inherit
// from qi::grammar, or they all contain a struct the inherits from grammar.  Right now
// VariableParser does the second, and Function Parser does the first.

#ifndef b2Test_grammar_h
#define b2Test_grammar_h

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/bind.hpp>

#include <memory>


#include <cmath>


#include "function_tree/node.hpp"
#include "function_tree/roots/function.hpp"

#include "function_tree/operators/arithmetic.hpp"
#include "function_tree/operators/trig.hpp"

#include "function_tree/symbols/number.hpp"
#include "function_tree/symbols/variable.hpp"




// see the following link for reading from arbitrary streams without putting the content into a std::string first
//http://boost-spirit.com/home/articles/qi-example/tracking-the-input-position-while-parsing/


// this solution for *lazy* make shared comes from the SO forum, user sehe.
// https://stackoverflow.com/questions/21516201/how-to-create-boost-phoenix-make-shared
//    post found using google search terms `phoenix construct shared_ptr`
namespace {
	template <typename T>
	struct make_shared_f
	{
		template <typename... A> struct result
		{ typedef std::shared_ptr<T> type; };

		template <typename... A>
		typename result<A...>::type operator()(A&&... a) const {
			return std::make_shared<T>(std::forward<A>(a)...);
		}
	};

	template <typename T>
	using make_shared_ = boost::phoenix::function<make_shared_f<T> >;
}







//////
//
//  the following code adapts the trig functions and other special functions to be phoenix-callable.
//
/////////

// http://www.boost.org/doc/libs/1_58_0/libs/phoenix/doc/html/phoenix/modules/function/adapting_functions.html

//BOOST_PHOENIX_ADAPT_FUNCTION(
//							 RETURN_TYPE
//							 , LAZY_FUNCTION
//							 , FUNCTION
//							 , FUNCTION_ARITY
//							 )



BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, cos_lazy, cos, 1);
BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, sin_lazy, sin, 1);
BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, tan_lazy, tan, 1);

BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, log_lazy, log, 1);
BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, exp_lazy, exp, 1);
BOOST_PHOENIX_ADAPT_FUNCTION(std::shared_ptr<bertini::node::Node>, sqrt_lazy, sqrt, 1);







namespace bertini {



	namespace qi = ::boost::spirit::qi;
	namespace ascii = ::boost::spirit::ascii;




	/**
	A Qi grammar parser for parsing text into function trees.  Currently called from the SystemParser.

	\todo Improve error detection and reporting for the FunctionParser.

	\brief A Qi grammar parser for parsing text into function trees.
	*/
	template<typename Iterator>
	struct FunctionParser : qi::grammar<Iterator, std::shared_ptr<node::Node>(), boost::spirit::ascii::space_type>
	{
		using Node = node::Node;
		using Function = node::Function;
		using Float = node::Float;
		using Integer = node::Integer;
		using Rational = node::Rational;
		
		FunctionParser(qi::symbols<char,std::shared_ptr<Node> > * encountered_symbols) : FunctionParser::base_type(root_rule_,"FunctionParser")
		{
			namespace phx = boost::phoenix;
			using qi::_1;
			using qi::_2;
			using qi::_3;
			using qi::_4;
			using qi::_val;
			using qi::eps;
			using qi::lit;

			using std::pow;
			using ::pow;

			root_rule_.name("function_");
			root_rule_ = expression_ [ _val = make_shared_<Function>()(_1)];


			///////////////////
			expression_.name("expression_");
			expression_ =
			term_ [_val = _1]
			>> *(   (lit('+') > term_ [_val += _1])
				 |  (lit('-') > term_ [_val -= _1])
				 )
			;

			term_.name("term_");
			term_ =
			factor_ [_val = _1]
			>> *(   (lit('*') > factor_ [_val *= _1])
				 |  (lit('/') > factor_ [_val /= _1])
				 )
			;

			factor_.name("factor_");
			factor_ =
			exp_elem_ [_val = _1]
			>> *(lit('^') // any number of ^somethings
				 > exp_elem_ [ phx::bind( []
										 (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
										 {
											 B = pow(B,P);
										 },
										 _val,_1)] )
			;

			exp_elem_.name("exp_elem_");
			exp_elem_ =
			(symbol_  >> !qi::alnum) [_val = _1]
			|   ( '(' > expression_  [_val = _1] > ')'  ) // using the > expectation here.
			|   (lit('-') > expression_  [_val = -_1])
			|   (lit('+') > expression_  [_val = _1])
			|   (lit("sin") > '(' > expression_ [_val = sin_lazy(_1)] > ')' )
			|   (lit("cos") > '(' > expression_ [_val = cos_lazy(_1)] > ')' )
			|   (lit("tan") > '(' > expression_ [_val = tan_lazy(_1)] > ')' )
			|   (lit("exp") > '(' > expression_ [_val = exp_lazy(_1)] > ')' )
			|   (lit("log") > '(' > expression_ [_val = log_lazy(_1)] > ')' )
			|   (lit("sqrt") > '(' > expression_ [_val = sqrt_lazy(_1)] > ')' )
			;








			symbol_.name("symbol_");
			symbol_ %=
			(*encountered_symbols) // the star here is the dereferencing of the encountered_symbols parameter to the constructor.
			|
			number_
			;











			number_.name("number_");
			number_ =
			long_number_string_ [ _val = make_shared_<Float>()(_1) ];



			long_number_string_.name("long_number_string_");
			long_number_string_ = eps[_val = std::string()] >>
			(
			 // 1. Read possible numbers before decimal, with possible negative
			 number_with_digits_after_point_ [_val += _1]
			 |
			 number_with_digits_before_point_ [_val += _1]
			 |
			 number_with_no_point_ [_val += _1]
			 >>   // reminder -- the - before the exponent_notation here means optional
			 - exponent_notation_ [_val+=_1]// Possible scientific notation, with possible negative in exponent.
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




			using qi::on_error;
			using boost::phoenix::val;
			using boost::phoenix::construct;


			on_error<qi::fail>
			(
			 root_rule_
			 , std::cout
			 << val("Function parser error:  expecting ")
			 << _4
			 << val(" here: \"")
			 << construct<std::string>(_3, _2)
			 << val("\"")
			 << std::endl
			 );





			//		debug(root_rule_);
			//		debug(expression_);
			//		debug(term_);
			//		debug(factor_);
			//		debug(exp_elem_);
			//		debug(number_);
			//		debug(number_with_no_point_);
			//		debug(number_with_digits_after_point_);
			//		debug(number_with_digits_before_point_);
			//		debug(exponent_notation_);
		}






		qi::rule<Iterator, std::shared_ptr<Node>(), ascii::space_type > root_rule_;
		// the rule for kicking the entire thing off

		qi::rule<Iterator, std::shared_ptr<Node>(), ascii::space_type> expression_, term_, factor_, exp_elem_;
		// rules for how to turn +-*/^ into operator nodes.



		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > symbol_;
		// any of the variables and numbers will be symbols.

		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > variable_;  // finds a previously encountered number, and associates the correct variable node with it.

		// the number_ rule wants to find strings from the various other number_ rules, and produces a Number node
		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > number_;

		// these rules all produce strings which are fed into numbers.
		qi::rule<Iterator, std::string()> long_number_string_, number_with_digits_before_point_, number_with_digits_after_point_, number_with_no_point_, exponent_notation_;
	};



} // re: namespace bertini

#endif
