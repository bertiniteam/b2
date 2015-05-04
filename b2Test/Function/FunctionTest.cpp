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
    num1->set<mpfr>(-3.3);
    num2->set<dbl>(4.3);
    num2->set<mpfr>(-4.3);
    
    std::unique_ptr<Node> add = std::make_unique<SumOperator>();
    dynamic_cast<SumOperator*>(add.get())->add_Child(std::move(num1), true);
    dynamic_cast<SumOperator*>(add.get())->add_Child(std::move(num2), false);

    std::unique_ptr<Node> num3 = std::make_unique<Constant>();
    std::unique_ptr<Node> num4 = std::make_unique<Constant>();
    num3->set<dbl>(3.3);
    num3->set<mpfr>(-3.3);
    num4->set<dbl>(4.3);
    num4->set<mpfr>(-4.3);

    std::unique_ptr<Node> mult= std::make_unique<MultOperator>();
    mult->add_Child(std::move(num3));
    mult->add_Child(std::move(num4));
    
//
    auto value = mult->eval<dbl>();
    for (int ii = 0; ii < 10; ++ii)
    {
        value = mult->eval<dbl>();
    }
//
    std::cout << value << std::endl;
    
    
    return 0;
}