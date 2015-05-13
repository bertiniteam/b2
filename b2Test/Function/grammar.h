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



#include "Node.h"
#include "sum_operator.h"
#include "mult_operator.h"
#include "negate_operator.h"
#include "exp_operator.h"
#include "constant.h"
#include "variable.h"


namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


// qi::grammar -> VariableParser
// This class describes the rule used to parse a list of variables and store a rule to
//  recognize variables in the function parse
template <typename Iterator>
class VariableParser
{
//    namespace phx = boost::phoenix;
public:
    // Constructor is used to define the grammar to parse variables.  Variables are separated by
    //  commas
    // TODO(JBC): We need to add support for underscore in variable names.
    
    
    
    // Accessor function for variable_
    qi::symbols<char,int> variable() { return variable_;};
    
    // Accessor for var_count
    int var_count() {return var_count_;};
    
    qi::rule<Iterator,boost::spirit::ascii::space_type> start() { return start_;};
   
    
private:
    // Method used to add variable names to the rule var_rule.  The index of the variable
    //  is stored as well, starting at 0.
    void AddVarToRule(std::vector<char> const& c)
    {
        std::string str;
        for(auto it : c)
        {
            str = str + it;
        }
        variable_.add(str,var_count_);
        var_count_++;
    }

    
    
    // Start rule used to parse variable list.
    qi::rule<Iterator,ascii::space_type> start_ = (*(qi::alnum | qi::char_('_')))[boost::phoenix::bind( [this](std::vector<char> const& c)
                                                                                        {
                                                                                            AddVarToRule(c);
                                                                                        }, qi::_1 ) ] % ",";

    
    // A rule defined to represent variables input by the user.  The rule recognizes
    //  whatever alphanumeric string(including underscore) used by the input file to define
    //  the variable, and returns the index of the variable, starting with 0.
    qi::symbols<char,int> variable_;
    
    
    
    // Keeps track of the index of the variable being parsed.
    int var_count_ = 0;
    
   
    
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
        sumexpr_ = (qi::eps)[_val = phx::new_<SumOperator>()] >>
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
                                          tempNeg->AddChild(std::shared_ptr<Node>(addinput));
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
                                          input->AddChild(std::shared_ptr<Node>(addinput));
                                      }
                                      ,_val, _1)] ) ) >>
        // 4. After done parsing, erase SumOp if there is only one term(i.e. no real summation happening)
        eps[ phx::bind( [](Node* &input)
                       {
                           int num = dynamic_cast<NaryOperator*>(input)->children_size();
                           if(num == 1)
                           {
                               Node* temp = dynamic_cast<NaryOperator*>(input)->first_child().get();
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
                               input->AddChild(std::shared_ptr<Node>(multinput));
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
                                          input->AddChild(std::shared_ptr<Node>(addinput));
                                      }
                                      ,_val, _1)] ) ) >>
        // 4. Erase MultOp if there is only one factor(i.e. no real multiplication happening)
        eps[ phx::bind( [](Node* &input)
                       {
                           int num = dynamic_cast<NaryOperator*>(input)->children_size();
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
                                      input->AddChild(std::shared_ptr<Node>(addinput));
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
        base_ =  mpfr_constant_[_val = phx::new_<Constant>(_1)]  |
        // 2. A variable possibly raise to a power with '^'.
        (
         // 2.a Before parsing, create _val = ExpOperator
         eps[_val = phx::new_<ExpOperator>()] >>
         // 2.b Add variable as base_ to _val
         variable_[ phx::bind( [](Node* input, int varnum, std::vector<std::shared_ptr<Variable> > input_vector)
                        {
                            input->AddChild(input_vector[varnum] );
                        }
                        ,_val, _1, *input_variables)] >>
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
         );
        /////////// NUMBER OR VARIABLE WITH EXP(base_) ////////////////

        
        
        
        
        
        /////////// MPFR Constant(mpfr_constant_) ////////////////
        mpfr_constant_ = eps[_val = std::string()] >>
        (
         // 1. Read possible numbers before decimal, with possible negative
         -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) >>
         // 2. Read possible numbers after the decimal
         -(qi::lit('.')[_val += "."]  >>  *(qi::char_(L'0',L'9')[_val += _1])) >>
         // 3. Possible scientific notation, with possible negative in exponent.
         -(qi::lit('e')[_val += "e"] >> -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) )
         );
        /////////// MPFR Constant(mpfr_constant_) ////////////////
    }
    
    
    
    
    
    
private:
    
    
    qi::rule<Iterator,Node*(),ascii::space_type> sumexpr_;
    qi::rule<Iterator,Node*(),ascii::space_type> multexpr_;
    qi::rule<Iterator,Node*(),ascii::space_type> subexpr_;
    qi::rule<Iterator,Node*(),ascii::space_type> parenexpr_;
    qi::rule<Iterator,Node*(),ascii::space_type> base_;
    qi::rule<Iterator,std::string()> mpfr_constant_;

};





#endif
