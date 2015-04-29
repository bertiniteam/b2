//
//  main.cpp
//  QiPlay
//
//  Created by Collins, James B. on 3/14/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
#define BOOST_SPIRIT_USE_PHOENIX_V3 1


#include <string>
#include <iostream>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/object/new.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>






namespace client
{
    
    
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    
    
    
    
    
    
    
    
    template<typename Iterator>
    struct MPFR : qi::grammar<Iterator,std::string()>
    {
        
        MPFR() : MPFR::base_type(MPnum)
        {
            namespace phx = boost::phoenix;
            using qi::_1;
            using qi::_val;
            using qi::eps;
            std::string s;
            
            MPnum = eps[_val = std::string()] >>
            (
             -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) >>
             -(qi::lit('.')[_val += "."]  >>  *(qi::char_(L'0',L'9')[_val += _1])) >>
             -(qi::lit('e')[_val += "e"] >> -(qi::lit('-')[_val += "-"]) >> *(qi::char_(L'0',L'9')[_val += _1]) )
            );
            
            //            sumexpr = eps[_val = phx::new_<Node>("add")] >>
            //            eps[phx::bind(&Node::addLeaf,_val,new Node("null",true))] >>
            //            var[phx::bind(&Node::print,_val)] >>
            //            *( ('-' >> var) |
            //              ('+' >> var) );
            
            
        }
        
        qi::rule<Iterator,std::string()> MPnum;
    };
    
    
    
}



int main()
{
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tMPFR Parser\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    
    typedef std::string::const_iterator iterator_type;
    typedef client::MPFR<iterator_type> MPnum;
    
    
    std::string str;
    
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your MPFR number ...or [q or Q] to quit\n\n";
    
    MPnum MP_parser;
    while (std::getline(std::cin, str))
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;
        
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        //[tutorial_roman_grammar_parse
        
        std::string test;
        bool r = phrase_parse(iter, end, MP_parser,client::ascii::space, test);
        if(test != "")
        {
            std::cout << "Decimal string is :" << test << std::endl;
        }
        else{
            std::cout << "didn't work\n";
        }
        
        
        if (r && iter == end)
        {
            std::cout << "-------------------------\n";
            std::cout << "Parsing succeeded\n";
            //            std::cout << "result = " << result << std::endl;
            std::cout << "-------------------------\n";
        }
        else
        {
            std::string rest(iter, end);
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << "stopped at: \" " << rest << "\"\n";
            std::cout << "-------------------------\n";
        }
        //]
    }
    
    std::cout << "Bye... :-) \n\n";
    return 0;
}


