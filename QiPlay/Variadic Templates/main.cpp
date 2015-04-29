//
//  main.cpp
//  Variadic Templates
//
//  Created by Collins, James B. on 4/23/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>

template<typename... TArgs>
class VariableStore
{
    
    
public:
    
    
    template<typename T>
    void set(T newVal)
    {
        get<T>() = newVal;
    }
    
    template<typename T>
    T& get()
    {
        return std::get<T>(vals);
    }
    
private:
    std::tuple<TArgs...> vals;
    
};





int main(int argc, const char * argv[]) {
    // insert code here...
    
    VariableStore<int, float, double> x;
    x.set<int>(4);
    x.set<float>(4.3f);
    x.set<double>(4.9);
    
    std::cout << x.get<int>() << std::endl
    << x.get<float>() << std::endl
    << x.get<double>() << std::endl;
    
    
    return 0;
}
