//
//  main.cpp
//  Variadic Templates
//
//  Created by Collins, James B. on 4/23/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>

template<typename DerivedType, typename... TArgs>
class Node
{
    
public:
    
    /* One eval function for each variable type */
    
    template<typename T>
    T eval()
    {
        return static_cast<DerivedType*>(this)->template eval<T>();
    }
    
    
    template<typename T>
    void set(T val)
    {
        std::get< std::pair<T,bool> >(this->current_value).first = val;
    }
    
    
//    virtual double eval(double) = 0;
//    virtual int eval(int) = 0;
    
    
    
    
protected:
    std::tuple< std::pair<TArgs,bool>... > current_value;
    
};



template<typename... TArgs>
class Sum : public Node<Sum<TArgs...>, TArgs...>
{
public:
//    Sum()
//    {
//        this->current_value = std::pair<float,bool> {7.9,true};
//    }
    
    
    
    
    template<typename T>
    T eval()
    {
        return std::get< std::pair<T,bool> >(this->current_value).first + 1;
    }
};




template<typename... TArgs>
class Mult : public Node<Mult<TArgs...>, TArgs...>
{
public:
//    Mult()
//    {
//        this->current_value = std::pair<float,bool> {7.9,true};
//    }
    
    
    template<typename T>
    T eval()
    {
        return std::get< std::pair<T,bool> >(this->current_value).first * 2;
    }
};




int main(int argc, const char * argv[]) {
    // insert code here...
    
    Node<Sum<float, double, int, long double>,float, double, int, long double>* x = new Sum<float, double, int, long double>();
    Node<Mult<float, double, int, long double>,float, double, int, long double>* y = new Mult<float, double, int, long double>();
    
    x->set<float>(1.2f);
    x->set<double>(3.1);
    x->set<int>(7);
    x->set<long double>(4.5676);
    
    y->set<float>(1.2f);
    y->set<double>(3.1);
    y->set<int>(7);
    y->set<long double>(4.5676);
    
    std::cout << x->eval<float>() << std::endl;
    std::cout << x->eval<double>() << std::endl;
    std::cout << x->eval<int>() << std::endl;
    std::cout << x->eval<long double>() << std::endl;
    
    std::cout << y->eval<float>() << std::endl;
    std::cout << y->eval<double>() << std::endl;
    std::cout << y->eval<int>() << std::endl;
    std::cout << y->eval<long double>() << std::endl;
    
    
    return 0;
}
