//
//  Grammar.h
//  b2Test
//
//  Created by Collins, James B. on 5/4/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_Grammar_h
#define b2Test_Grammar_h

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




namespace parser
{
    
    
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    
    
    
    struct var_ : qi::symbols<char,int>
    {
        
    }var;
    
    
    
    int varcount = 0;
    void addVars(std::vector<char> const& c)
    {
        std::string str;
        for(auto it : c)
        {
            str = str + it;
        }
        var.add(str,varcount);
        varcount++;
    }
    
    template<typename Iterator>
    struct variables : qi::grammar<Iterator,ascii::space_type>
    {
        
        variables() : variables::base_type(start)
        {
            namespace phx = boost::phoenix;
            using qi::_1;
            
            start = (*qi::alnum)[&addVars] % ",";
        }
        
        qi::rule<Iterator,ascii::space_type> start;
    };

    template<typename Iterator>
    struct polynomial : qi::grammar<Iterator,Node*(),ascii::space_type>
    {
        
        qi::rule<Iterator,Node*(),ascii::space_type> sumexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> multexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> subexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> parenexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> base;

        
        polynomial() : polynomial::base_type(sumexpr)
        {
            namespace phx = boost::phoenix;
            using qi::_1;
            using qi::_val;
            using qi::eps;
            
            
            /////////// TERMS(sumexpr) ////////////////
            // Parses a list of terms(see FACTORS for def) separated by '+' or '-'
            
            // 1. Before start parsing, create a SumOperator
            sumexpr = (qi::eps)[_val = phx::new_<SumOperator>()] >>
            // 2.a Add first term to _val
            ( multexpr[ phx::bind( [](Node* input, Node* addinput)
                           {
                               dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), true);
                           }
                           ,_val, _1)] |
             // 2.b Create Negate and add first term to it, add Negate to _val
             ( '-'>> multexpr)[ phx::bind( [](Node* input, Node* addinput)
                                          {
                                              std::shared_ptr<NegateOperator> tempNeg = std::make_shared<NegateOperator>();
                                              tempNeg->AddChild(std::shared_ptr<Node>(addinput));
                                              dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(tempNeg));
                                          }
                                          ,_val, _1)]) >>
            // 3. Add other terms to _val
            // 3.a If '-' in front, add as negative term
            *( ('-' >> multexpr[ phx::bind( [](Node* input, Node* addinput)
                                           {
                                               dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), false);
                                           }
                                           ,_val, _1)] ) |
            // 3.b If '+' in front, add as positive term
              ('+' >> multexpr[ phx::bind( [](Node* input, Node* addinput)
                                          {
                                              dynamic_cast<SumOperator*>(input)->AddChild(std::shared_ptr<Node>(addinput), true);
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

            
            
            
            
            /////////// FACTORS(multexpr) ////////////////
            // Parses a list of factors(see BASE OR PARENS for def) separated by '*'
            
            // 1. Before parsing, create a _val = MultOperator
            multexpr = eps[_val = phx::new_<MultOperator>()] >>
            // 2. Add first factor to _val
            subexpr[ phx::bind( [](Node* input, Node* multinput)
                               {
                                   input->AddChild(std::shared_ptr<Node>(multinput));
                               }
                               ,_val, _1)] >>
            // 3. Parse each factor and add to _val
            *( '*' >> subexpr )[ phx::bind( [](Node* input, Node* multinput)
                                           {
                                               input->AddChild(std::shared_ptr<Node>(multinput));
                                           }
                                           ,_val, _1)] >>
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

            
            /////////// BASE OR PARENS(subexpr) ////////////////
            // Either a symbol(see base for def) or something in parenthesis.
            // NOTE: Can be anything in parenthesis
            subexpr = base[_val = _1] | parenexpr[_val = _1];
            
            
            /////////// PARENS WITH EXP ////////////////
            // Anything in parenthesis possibly raised to power with '^'.  Uses recurence with sumexpr to parse
            // expression within parenthesis.
            
            // 1. Before parsing, create _val = ExpOperator
            parenexpr = eps[_val = phx::new_<ExpOperator>()] >>
            // 2. Parse something inside parens and add to _val
            '(' >> sumexpr[ phx::bind( [](Node* input, Node* addinput)
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
            
            
            
            
            
            /////////// NUMBER OR VARIABLE WITH EXP(base) ////////////////
            // These are the symbols or leaves in the function tree.
            
            // 1. Any double number.  TODO(JBC): Change with MPFR parser!!!
            base =  qi::double_[_val = phx::new_<Constant>(_1)]  |
            // 2. A variable possibly raise to a power with '^'.
            (
             // 2.a Before parsing, create _val = ExpOperator
             eps[_val = phx::new_<ExpOperator>()] >>
             // 2.b Add variable as base to _val
             var[ phx::bind( [](Node* input, int varnum)
                            {
                                input->AddChild(std::make_shared<Variable>(varnum));
                            }
                            ,_val, _1)] >>
             // 2.c Possibly parse exponent with '^' notation and set exp in _val
             -('^' >> qi::int_[ phx::bind( [](Node* input, int expinput)
                                          {
                                              dynamic_cast<ExpOperator*>(input)->set_exponent(expinput);
                                          }
                                          ,_val, _1)]) >>
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
        }
    };
}
#endif
