//
//  FunctionTest.cpp
//  b2Test
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>

#include "Node.h"
#include "SumOperator.h"
#include "MultOperator.h"
#include "Constant.h"



int main()
{
    std::cout << "Hello World!\n";
    
    std::unique_ptr<Node> num1 = std::make_unique<Constant>();
    std::unique_ptr<Node> num2 = std::make_unique<Constant>();
    num1->set<dbl>(3.3);
    num2->set<dbl>(2.3);
    
    std::unique_ptr<Node> add = std::make_unique<SumOperator>();
    add->add_Child(std::move(num1));
    add->add_Child(std::move(num2));
//
    auto value = add->eval<dbl>();
//
    std::cout << value << std::endl;
    
    
    return 0;
}