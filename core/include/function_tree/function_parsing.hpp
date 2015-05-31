//This file is part of Bertini 2.0.
//
//Foobar is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Foobar is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// grammar.h:  This file contains all the grammars used to parse an input file.
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


#include "function_tree/Node.h"
#include "function_tree/roots/function.h"
#include "function_tree/operators/sum_operator.h"
#include "function_tree/operators/mult_operator.h"
#include "function_tree/operators/negate_operator.h"
#include "function_tree/operators/exp_operator.h"
#include "function_tree/symbols/number.h"
#include "function_tree/symbols/variable.h"

#include "function_tree/operators/binary_operator.h"



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
			valid_variable_name_ %= +qi::alpha >> *(qi::alnum | qi::char_('_'));
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
	struct BrakeParser : qi::grammar<Iterator, std::shared_ptr<Node>(), boost::spirit::ascii::space_type>
	{
		
		//,"BRAKEPARSER"
		// sumexpr is going to be the first rule called, the start rule.
		BrakeParser(const std::vector<std::shared_ptr<Variable> > * input_variables, qi::symbols<char,int> encountered_variables) : BrakeParser::base_type(root_rule_,"FunctionParser")
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
			symbol_ [_val = _1]
			|   ( '(' > expression_  [_val = _1] > ')'  ) // using the > expectation here.
			|   (lit('-') > expression_  [_val = -_1])
			|   (lit('+') > expression_  [_val = _1])
			;
			
			
			
			
			
			
			
			
			symbol_.name("symbol_");
			symbol_ %=
			variable_
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
		
		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > variable_;  // finds a previously encountered number, and associates the correct variable node with it.
		
		// the number_ rule wants to fund strings from the various other number_ rules, and produces a Number node
		qi::rule<Iterator, std::shared_ptr<Node>(),  ascii::space_type > number_;
		
		// these rules all produce strings which are fed into numbers.
		qi::rule<Iterator, std::string()> long_number_string_, number_with_digits_before_point_, number_with_digits_after_point_, number_with_no_point_, exponent_notation_;
	};
	
	
	
	
	
	//qi::grammar -> FunctionParser
	// This class describes the grammar used to parse a function.  The variable list must be passed
	//  to the constructor.  During parsing, the function tree is built by the grammar.
	template <typename Iterator>
	struct FunctionParser : qi::grammar<Iterator,Node*(),boost::spirit::ascii::space_type>
	{
		
		
	public:
		FunctionParser(const std::vector<std::shared_ptr<Variable> >* input_variables, qi::symbols<char,int> variable_) :
		FunctionParser::base_type(sumexpr_)
		{
			namespace phx = boost::phoenix;
			using qi::_1;
			using qi::_val;
			using qi::eps;
			
			
			
			
			//      TODO:(JBC) Currently the function parser is a factory for creating raw pointers using phoenix(phx::new_).
			// We need to check is we can use shared_ptr with Qi.  If not, check if we can use C++ new as we know more of how
			// that works.
			//      If tried having Qi return shared_ptr<Node> and it wouldn't compile.  I may have missed something though.
			// I've also tried using C++ new and Qi tried to create two copies of everything.  That was before I used lambdas
			// to do actions though.  With lambda, the C++ new might work better
			//
			//      TODO:(JBC) If SumOperator has only one term, I delete the SumOperator node and just keep the term.  I do
			// this with a raw delete, and I'm pretty sure that's not good.  Need to look into this and see if there is a
			// better way.
			//
			//      TODO:(JBC) When everything is up and running, need to check if compile times are prohibitively long due to
			// the many lambdas and actions done within the parser itself.  Some options that can be taken out are deleting
			// a SumOperator or MultOperator when only a single term or factor is present.  This could be done in post-processing if
			// compile times are getting too long.
			//
			//      TODO:(JBC) Generalize the parser to be able to read anything Matlab can read.  For example, we need to allow
			// the parser to read 4 +-9 as -5.  Right now multiple operators together cannot be parsed.  Determine all other
			// possibilities allowed by Matlab and add them to the parser.
			
			
			/////////// TERMS(sumexpr_) ////////////////
			// Parses a list of terms(see FACTORS for def) separated by '+' or '-'
			
			// 1. Before start parsing, create a SumOperator
			sumexpr_ = qi::eps[_val = phx::new_<SumOperator>()] >>
			// 2.a Add first term to _val
			( multexpr_[ phx::bind( [](Node* input, Node* addinput)
								   {
									   dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), true);
								   }
								   ,_val, _1)] |
			 // 2.b Create Negate and add first term to it, add Negate to _val
			 ( '-'>> multexpr_)[ phx::bind( [](Node* input, Node* addinput)
										   {
											   std::shared_ptr<NegateOperator> tempNeg = std::make_shared<NegateOperator>();
											   tempNeg->SetChild(std::shared_ptr<Node>(addinput));
											   dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(tempNeg));
										   }
										   ,_val, _1)] ) >>
			// 3. Add other terms to _val
			// 3.a If '-' in front, add as negative term
			*( ('-' >> multexpr_[ phx::bind( [](Node* input, Node* addinput)
											{
												dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), false);
											}
											,_val, _1)] ) |
			  // 3.b If '+' in front, add as positive term
			  ('+' >> multexpr_[ phx::bind( [](Node* input, Node* addinput)
										   {
											   dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput));
										   }
										   ,_val, _1)] ) ) >>
			// 4. After done parsing, erase SumOp if there is only one term(i.e. no real summation happening)
			eps[ phx::bind( [](Node* &input)
						   {
							   auto num = dynamic_cast<NaryOperator*>(input)->children_size();
							   if(num == 1)
							   {
								   Node* temp = dynamic_cast<NaryOperator* >(input)->first_child().get();
								   std::cout << "A" << std::endl;
								   delete input;  //This is dangerous!!!!!!
								   input = temp;
							   }
						   },_val)];
			/////////// TERMS(sumexpr_) ////////////////
			
			
			
			
			/////////// FACTORS(multexpr_) ////////////////
			// Parses a list of factors(see BASE OR PARENS for def) separated by '*'
			
			// 1. Before parsing, create a _val = MultOperator
			multexpr_ = eps[_val = phx::new_<MultOperator>()] >>
			// 2. Add first factor to _val
			subexpr_[ phx::bind( [](Node* input, Node* multinput)
								{
									dynamic_cast<MultOperator*>(input)->AddChild(std::shared_ptr<Node>(multinput));
								}
								,_val, _1)] >>
			// 3. Add other factors to _val
			// 3.a If '/' in front, divide factor
			*( ('/' >> subexpr_[ phx::bind( [](Node* input, Node* addinput)
										   {
											   dynamic_cast<MultOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), false);
										   }
										   ,_val, _1)] ) |
			  // 3.b If '*' in front, add as positive term
			  ('*' >> subexpr_[ phx::bind( [](Node* input, Node* addinput)
										  {
											  dynamic_cast<MultOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput));
										  }
										  ,_val, _1)] ) ) >>
			// 4. Erase MultOp if there is only one factor(i.e. no real multiplication happening)
			eps[ phx::bind( [](Node* &input)
						   {
							   auto num = dynamic_cast<NaryOperator*>(input)->children_size();
							   if(num == 1)
							   {
								   Node* temp = dynamic_cast<NaryOperator*>(input)->first_child().get();
								   delete input;  //This is dangerous!!!!!!
								   input = temp;
							   }
						   },_val)];
			/////////// FACTORS(multexpr_) ////////////////
			
			
			
			
			
			
			/////////// BASE OR PARENS(subexpr_) ////////////////
			// Either a symbol(see base_ for def) or something in parenthesis.
			// NOTE: Can be anything in parenthesis
			subexpr_ = base_[_val = _1] | parenexpr_[_val = _1];
			/////////// BASE OR PARENS(subexpr_) ////////////////
			
			
			
			
			
			
			/////////// PARENS WITH EXP(parenexpr_) ////////////////
			// Anything in parenthesis possibly raised to power with '^'.  Uses recurence with sumexpr_ to parse
			// expression within parenthesis.
			
			// 1. Before parsing, create _val = ExpOperator
			parenexpr_ = eps[_val = phx::new_<ExpOperator>()] >>
			// 2. Parse something inside parens and add to _val
			'(' >> sumexpr_[ phx::bind( [](Node* input, Node* addinput)
									   {
										   dynamic_cast<ExpOperator*>(input)->SetChild(std::shared_ptr<Node>(addinput));
									   }
									   ,_val, _1)] >> ')' >>
			// 3. Parse the exponent and set the exponent in _val.  If no exp to parse
			//     exp in Operator is set to 1 by default
			-( '^' >> qi::int_[ phx::bind( [](Node* input, int expinput)
										  {
											  dynamic_cast<ExpOperator*>(input)->set_exponent(expinput);
										  }
										  ,_val, _1)] ) >>
			// 4. Erase ExpOp if exp = 1(i.e. no real exponentiation happening)
			eps[ phx::bind( [](Node* &input)
						   {
							   int exp = dynamic_cast<ExpOperator*>(input)->exponent();
							   if(exp == 1)
							   {
								   Node* temp = dynamic_cast<UnaryOperator*>(input)->first_child().get();
								   delete input;  //This is dangerous!!!!!!
								   input = temp;
							   }
						   },_val)];
			/////////// PARENS WITH EXP(parenexpr_) ////////////////
			
			
			
			
			
			
			
			/////////// NUMBER OR VARIABLE WITH EXP(base_) ////////////////
			// These are the symbols or leaves in the function tree.
			
			// 1. Any double number.  TODO(JBC): Change with MPFR parser!!!
			base_ =
			// 2. A variable possibly raise to a power with '^'.
			(
			 // 2.a Before parsing, create _val = ExpOperator
			 eps[_val = phx::new_<ExpOperator>()] >>
			 // 2.b Add variable as base_ to _val
			 variable_[ phx::bind( [](Node* input, int varnum, const std::vector<std::shared_ptr<Variable> > * input_vector)
								  {
									  dynamic_cast<ExpOperator*>(input)->SetChild((*input_vector)[varnum] );
								  }
								  ,_val, _1, input_variables)] >>
			 // 2.c Possibly parse exponent with '^' notation and set exp in _val
			 -( '^' >> qi::int_[ phx::bind( [](Node* input, int expinput)
										   {
											   dynamic_cast<ExpOperator*>(input)->set_exponent(expinput);
										   }
										   ,_val, _1)] ) >>
			 // 2.d Erase ExpOp if exp = 1(i.e. no real exponentiation happening)
			 eps[ phx::bind( [](Node* &input)
							{
								int exp = dynamic_cast<ExpOperator*>(input)->exponent();
								if(exp == 1)
								{
									Node* temp = dynamic_cast<UnaryOperator*>(input)->first_child().get();
									delete input;  //This is dangerous!!!!!!
									input = temp;
								}
							},_val)]
			 ) |
			mpfr_constant_[_val = phx::new_<Number>(_1)]  ;
			/////////// NUMBER OR VARIABLE WITH EXP(base_) ////////////////
			
			
			
			
			
			
			/////////// MPFR Number(mpfr_constant_) ////////////////
			mpfr_constant_ = eps[_val = std::string()] >>
			(
			 // 1. Read possible numbers before decimal, with possible negative
			 -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) >>
			 // 2. Read possible numbers after the decimal
			 -(qi::lit('.')[_val += "."]  >>  *(qi::char_(L'0',L'9')[_val += _1])) >>
			 // 3. Possible scientific notation, with possible negative in exponent.
			 -(qi::lit('e')[_val += "e"] >> -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) )
			 );
			/////////// MPFR Number(mpfr_constant_) ////////////////
		}
		
		
		
		
		
		
	private:
		
		
		qi::rule<Iterator,Node*(),ascii::space_type> sumexpr_;
		qi::rule<Iterator,Node*(),ascii::space_type> multexpr_;
		qi::rule<Iterator,Node*(),ascii::space_type> subexpr_;
		qi::rule<Iterator,Node*(),ascii::space_type> parenexpr_;
		qi::rule<Iterator,Node*(),ascii::space_type> base_;
		qi::rule<Iterator,std::string()> mpfr_constant_;
		
	};
	
	
	
} // re: namespace bertini

#endif
