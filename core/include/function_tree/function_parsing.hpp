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


#include "function_tree/Node.hpp"
#include "function_tree/roots/function.hpp"
#include "function_tree/operators/sum_operator.hpp"
#include "function_tree/operators/mult_operator.hpp"
#include "function_tree/operators/negate_operator.hpp"
#include "function_tree/operators/expon_operator.hpp"
#include "function_tree/symbols/number.hpp"
#include "function_tree/symbols/variable.hpp"

#include "function_tree/operators/binary_operator.hpp"



// see the following link for reading from arbitrary streams without putting the content into a std::string first
//http://boost-spirit.com/home/articles/qi-example/tracking-the-input-position-while-parsing/


// this solution for lazy make shared comes from the SO forum, user sehe.
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












// the following code is adapted from
// http://boost.2283326.n4.nabble.com/shared-ptr-how-to-interact-with-shared-ptr-td4665327.html

//namespace boost { namespace spirit { namespace traits
//	{
//		template<class T> using shared_ptr = std::shared_ptr<T>;
//		using std::make_shared;
//
//
//
//		namespace detail {
//			template <typename T>
//			struct transform_to_sharedptr : boost::static_visitor<shared_ptr<T> > {
//				shared_ptr<T> operator()(shared_ptr<T>  const& v) const { return v; }
//
//				template <typename U>
//				shared_ptr<T> operator()(shared_ptr<U>  const& v) const { return v; }
//
//				template <typename... Ts>
//				shared_ptr<T> operator()(boost::variant<Ts...> const& v) const { return boost::apply_visitor(*this, v); }
//
//				template <typename U>
//				shared_ptr<T> operator()(U                     const& v) const { return make_shared<U>(v); }
//			};
//		}
//
//		template<typename T, typename U>
//		struct assign_to_attribute_from_value<shared_ptr<T>, U> {
//			static void call(U const& v, shared_ptr<T>& attr) {
//				std::cout << __PRETTY_FUNCTION__ << "\n";
//				attr = detail::transform_to_sharedptr<T>()(v);
//			}
//		};
//	}}}




namespace bertini {
	
	
	
	namespace qi = ::boost::spirit::qi;
	namespace ascii = ::boost::spirit::ascii;
	
	
	// qi::grammar -> VariableParser
	// This class describes the rule used to parse a list of variables and store a rule to
	//  recognize variables in the function parse
	template <typename Iterator>
	class VariableParser : qi::grammar<Iterator, boost::spirit::ascii::space_type>
	{
		//    namespace phx = boost::phoenix;
	public:
		
		// Constructor is used to define the grammar to parse variables.  Variables are separated by
		//  commas
		VariableParser() : VariableParser::base_type(start_,"VariableParser")
		{
			start_.name("VariableParser.start_");
			start_ =
			valid_variable_name_ [boost::phoenix::bind( [this](std::string const& c)
															 {
																 AddVarToRule(c);
															 }, qi::_1 ) ]
			% ","; // variables are separated by commas
			
			valid_variable_name_.name("valid_variable_name_");
            valid_variable_name_ = +qi::alpha >> *(qi::alnum | qi::char_('_')) >> -('[' >> +qi::digit >> ']')  ;
            
		}
		
		
		
		
		qi::rule<Iterator,boost::spirit::ascii::space_type> start() { return start_;};
		
		
		
		
		// Accessor function for variable_
		qi::symbols<char,int> variable() { return variable_;};
		
		// Accessor for var_count
		int var_count() {return var_count_;};
		
		
		
		std::vector<std::shared_ptr<Variable> > get_var_group(){
			return produced_variables_;
		}
		
	private:
		
		
		// Method used to add variable names to the rule var_rule.  The index of the variable
		//  is stored as well, starting at 0.
		void AddVarToRule(std::string const& c)
		{
			produced_variables_.push_back(std::make_shared<Variable>(c));
			
			variable_.add(c,var_count_);
			var_count_++;
		}
		
		
		
		
		// Start rule used to parse variable list.
		qi::rule<Iterator,ascii::space_type> start_;
		
		// the rule which determines valid variable names
		qi::rule<Iterator, std::string()> valid_variable_name_;
        
		
		
		
		// A rule defined to represent variables input by the user.  The rule recognizes
		//  whatever alphanumeric string(including underscore) used by the input file to define
		//  the variable, and returns the index of the variable, starting with 0.
		qi::symbols<char,int> variable_;
		
		
		// vector of shared pointers to nodes which represent a variable group or whatever.
		std::vector< std::shared_ptr<Variable> > produced_variables_;
		
		
		// Keeps track of the index of the variable being parsed.
		int var_count_ = 0;
		
	};
	
	
	
	
	
	
	
	template<typename Iterator>
	struct FunctionParser : qi::grammar<Iterator, std::shared_ptr<Node>(), boost::spirit::ascii::space_type>
	{
		
		//,"BRAKEPARSER"
		// sumexpr is going to be the first rule called, the start rule.
		FunctionParser(const std::vector<std::shared_ptr<Variable> > * input_variables, qi::symbols<char,int> encountered_variables) : FunctionParser::base_type(root_rule_,"FunctionParser")
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
			
			
			root_rule_.name("root_rule_");
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
            ( '(' > expression_  [_val = _1] > ')'  ) // using the > expectation here.
			|   (lit('-') > expression_  [_val = -_1])
			|   (lit('+') > expression_  [_val = _1])
            |   (lit("sin(") > expression_ [ phx::bind( []
                                                       (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
                                                       {
                                                           B = sin(P);
                                                       },
                                                       _val,_1)] > lit(')') )
            |   (lit("cos(") > expression_ [ phx::bind( []
                                                       (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
                                                       {
                                                           B = cos(P);
                                                       },
                                                       _val,_1)] > lit(')') )
            |   (lit("tan(") > expression_ [ phx::bind( []
                                                       (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
                                                       {
                                                           B = tan(P);
                                                       },
                                                       _val,_1)] > lit(')') )
            |   (lit("exp(") > expression_ [ phx::bind( []
                                                       (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
                                                       {
                                                           B = exp(P);
                                                       },
                                                       _val,_1)] > lit(')') )
            |   (lit("sqrt(") > expression_ [ phx::bind( []
                                                       (std::shared_ptr<Node> & B, std::shared_ptr<Node> P)
                                                       {
                                                           B = sqrt(P);
                                                       },
                                                       _val,_1)] > lit(')') )
            |   symbol_ [_val = _1]
			;
			
			
			
			
			
			
			
			
			symbol_.name("symbol_");
			symbol_ %=
			variable_
			|
			constant_
            |
            number_
			// the following commented out lines are place holders for incorporation later.
			//		|
			//		constant_
			//		|
			//		subfunction_
			//		|
			//		parameter_
			;
			
			
			variable_.name("variable_");
			variable_ =
			encountered_variables [ phx::bind( [](std::shared_ptr<Node> & input, int varnum, const std::vector<std::shared_ptr<Variable> > * input_vector)
											  {
												  input = (*input_vector)[varnum];
											  }
											  ,_val, _1, input_variables)];
			
			
			
			
            constant_.name("constant_");
            constant_ =
            qi::lit("pi") [_val = make_shared_<SpecialNumber>()("pi")]
            |
            qi::lit("Pi") [_val = make_shared_<SpecialNumber>()("pi")]
            |
            qi::lit('e') [_val = make_shared_<SpecialNumber>()("e")]
            |
            qi::lit('E') [_val = make_shared_<SpecialNumber>()("e")]
            |
            qi::lit('i') [_val = make_shared_<SpecialNumber>()("i")]
            |
            qi::lit('I') [_val = make_shared_<SpecialNumber>()("i")]
            |
            qi::lit("1i") [_val = make_shared_<SpecialNumber>()("i")];
            
            
			number_.name("number_");
			number_ =
			long_number_string_ [ _val = make_shared_<Number>()(_1) ];
			
			
			
			long_number_string_.name("long_number_string_");
			long_number_string_ = eps[_val = std::string()] >>
			(
			 // 1. Read possible numbers before decimal, with possible negative
			 number_with_digits_after_point_ [_val += _1]
			 |
			 number_with_digits_before_point_ [_val += _1]
			 |
			 number_with_no_point_ [_val += _1]
			 >   // reminder -- the - before the exponent_notation here means optional
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
			-(qi::char_(L'0',L'9')[_val += _1]) // find any number of digits after the point
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
			 << val("Error! Expecting ")
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
		
		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > variable_;  // finds a previously encountered variable, and associates the correct variable node with it.
		
		// the number_ rule wants to fund strings from the various other number_ rules, and produces a Number node
		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > number_;

        // the constant_ rule wants to find strings representing a special constant, and produces a SpecialNumber node
        qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > constant_;

		// these rules all produce strings which are fed into numbers.
		qi::rule<Iterator, std::string()> long_number_string_, number_with_digits_before_point_, number_with_digits_after_point_, number_with_no_point_, exponent_notation_;
	};
	
	
	
	
	
	
} // re: namespace bertini

#endif
