//
//  FunctionTest.cpp
//  b2Test
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>
#include <vector>

#include "variable.h"


// Obviously a bad idea :-)
std::vector<std::shared_ptr<Variable> > variable_nodes;



#include "grammar.h"

#include "node.h"
#include "sum_operator.h"
#include "mult_operator.h"
#include "negate_operator.h"
#include "constant.h"





////////// TESTING /////////////
int Node::tabcount = 0;
////////// TESTING /////////////


int main()
{
    std::cout.precision(10);
    
    // 1. Read in the variables //
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tPolynomial Parser\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your symbols separated by commas ...or [q or Q] to quit\n\n";
    
    typedef std::string::const_iterator iterator_type;
    typedef parser::variables<iterator_type> variables;
    typedef parser::polynomial<iterator_type> polys;
    
    variables var_parser;
    
    std::string str;
    std::getline(std::cin, str);
    if (str.empty() || str[0] == 'q' || str[0] == 'Q')
        return 0;
    
    iterator_type it = str.begin();
    iterator_type end = str.end();
    phrase_parse(it,end,var_parser,parser::ascii::space);
    
    for (int ii = 0; ii < parser::varcount; ++ii)
    {
        double variable_input;
        std::cout << "Variable " << ii << " = ";
        std::getline(std::cin,str);
        variable_input = std::stod(str);
        Variable new_variable(variable_input);
        variable_nodes.push_back(std::make_shared<Variable>(variable_input) );
    }
    
    
    
    // 2. Main loop //
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your polynomial ...or\n";
    std::cout << "Enter e to evaluate a polynomial...or\n";
    std::cout << "Enter s to set the variables again...or\n";
    std::cout << "Enter [q or Q] to quit\n\n";
    
    polys poly_parser;
    Node* test_function{nullptr};
    while (std::getline(std::cin, str))
    {
        if (str[0] == 'q' || str[0] == 'Q')
        {
            break;
        }
        else if(str[0] == 'e')
        {
            if(test_function != nullptr)
            {
                std::cout << "test_function  =  " << test_function->Eval<dbl>() << std::endl;
                std::cout << "test_function tree is:\n";
                test_function->PrintTree();
            }
            else
            {
                std::cout << "No function defined!\n";
            }
            continue;
        }
        else if(str[0] == 's')
        {
            for (int ii = 0; ii < parser::varcount; ++ii)
            {
                double variable_input;
                std::cout << "Variable " << ii << " = ";
                std::getline(std::cin,str);
                variable_input = std::stod(str);
                variable_nodes[ii]->set_current_value<dbl>(variable_input);
            }
            test_function->Reset<dbl>();
            continue;
        }
        
        
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();

        bool r = phrase_parse(iter, end, poly_parser,parser::ascii::space, test_function);
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
        
        
        std::cout << "test_function  =  " << test_function->Eval<dbl>() << std::endl;
        std::cout << "test_function tree is:\n";
        test_function->PrintTree();
    }
    
    return 0;
}