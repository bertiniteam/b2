//
//  Node.h
//  b2Test
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//




// Node class
// Descrip: Abtract base class for all operators and children in a function tree


#ifndef __b2Test__Node__
#define __b2Test__Node__




#include <iostream>
#include <string>
#include <tuple>

#include <boost/type_index.hpp>

using dbl = double;
using mpfr = long double;



class Node
{
protected:
    std::tuple< std::pair<dbl,bool>, std::pair<mpfr,bool> > current_value;
    
    
    
    ///////// OVERRIDEN METHODS /////////////////
    virtual dbl fresh_eval(dbl) = 0;
    virtual mpfr fresh_eval(mpfr) = 0;
    
    ///////// OVERRIDEN METHODS /////////////////
    
    
public:
    
    
    ////////////TESTING////////////////    
    static int tabcount;
    virtual void printTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << boost::typeindex::type_id_runtime(*this).pretty_name() << std::endl;
    }

    ////////////TESTING////////////////
    
    
    ///////// OVERRIDEN METHODS /////////////////
    virtual std::string get_string() = 0;
    virtual void add_Child(std::shared_ptr<Node> child) = 0;
    virtual void add_ChildQi(Node* child) {}; // We need this b/c unique_ptrs don't work well with Qi

    ///////// OVERRIDEN METHODS /////////////////
    
    template <typename T>
    void set(T value)
    {
        std::get< std::pair<T,bool> >(current_value).first = value;
        std::get< std::pair<T,bool> >(current_value).second = false;
    }
    
    
    
    
    //Evaluate the node.  If flag false, just return value, if flag true
    //run the specific fresh_eval of the node
    template<typename T>
    T eval()
    {
        auto& val_pair = std::get< std::pair<T,bool> >(current_value);
        if(!val_pair.second)
        {
//            std::cout << "Fresh eval\n";
            T input{};
            val_pair.first = fresh_eval(input);
            val_pair.second = true;
            return val_pair.first;
        }
        else
        {
            return val_pair.first;
        }
    }
    

    
    //Should we have a method to reset all values to fresh eval?
    
    //reset only double values to fresh eval
    template <typename T>
    void reset()
    {
        std::get< std::pair<T,bool> >(current_value).second = false;
    }
    
    

};



#endif /* defined(__b2Test__Node__) */
