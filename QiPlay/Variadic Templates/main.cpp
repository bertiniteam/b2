//
//  main.cpp
//  Variadic Templates
//
//  Created by Collins, James B. on 4/23/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>

template<typename... TArgs>
struct VariableStore
{
    std::tuple<TArgs...> vals;
    
    template<typename T>
    void add(T newVal)
    {
        get<T>().push_back(newVal);
    }
    
    template<typename T>
    T& get()
    {
        return std::get<T>(vals);
    }
    
};





int main(int argc, const char * argv[]) {
    // insert code here...
    
    VariableStore<int, float, double> x;
    x.add<int>(4);
    x.add<float>(4.3f);
    x.add<double>(4.9);
    
    std::cout << x.get<int>() << std::endl
    << x.get<float>() << std::endl
    << x.get<double>() << std::endl;
    
    
    return 0;
}
