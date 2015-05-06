//
//  ParseTest.cpp
//  b2Test
//
//  Created by Collins, James B. on 5/4/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>
#include <string>

#include "Grammar.h"
#include "Node.h"
#include "SumOperator.h"



int main()
{
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tPolynomial Parser\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your symbols separated by commas ...or [q or Q] to quit\n\n";
    
    typedef std::string::const_iterator iterator_type;
    typedef parser::variables<iterator_type> variables;
    typedef parser::polynomial<iterator_type> polys;
    
    variables var_parser; // Our grammar
    
    std::string str;
    std::getline(std::cin, str);
    if (str.empty() || str[0] == 'q' || str[0] == 'Q')
        return 0;
    
    iterator_type it = str.begin();
    iterator_type end = str.end();
    phrase_parse(it,end,var_parser,parser::ascii::space);
    
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your polynomial ...or [q or Q] to quit\n\n";
    
    polys poly_parser;
    Node* test{nullptr};
    while (std::getline(std::cin, str))
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;
        
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        //[tutorial_roman_grammar_parse
        
        bool r = phrase_parse(iter, end, poly_parser,parser::ascii::space, test);
       
        
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
            std::cout << "stopped at: \": " << rest << "\"\n";
            std::cout << "-------------------------\n";
        }
        //]
    }
    
    
    if(test != nullptr)
    {
        std::unique_ptr<Node> test2(test);
    }
    
    std::cout << "Bye... :-) \n\n";
    return 0;
}
