//
//  main.cpp
//  Variadic Templates
//
//  Created by Collins, James B. on 4/23/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>

template<typename... TArgs>
class Node
{
    
public:
    
    /* One eval function for each variable type */
    float eval(float){return 2.0;};
//    virtual double eval(double) = 0;
//    virtual int eval(int) = 0;
    
    
    
    
protected:
    std::tuple< std::pair<TArgs,bool>... > current_value;
    
};



template<typename... TArgs>
class Sum : public Node<TArgs...>
{
public:
    Sum()
    {
        this->current_value = std::pair<float,bool> {7.4,true};
    }
    
    
    float eval(float)
    {
        return std::get< std::pair<float,bool> >(this->current_value).first;
    }
};





int main(int argc, const char * argv[]) {
    // insert code here...
    
    Node<float>* x = new Sum<float>();
    
    std::cout << x->eval(8.7f) << std::endl;
    
    
    return 0;
}
