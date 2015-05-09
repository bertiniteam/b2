// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// node.h:  Declares the class Node.


#ifndef __b2Test__Node__
#define __b2Test__Node__




#include <iostream>
#include <string>
#include <tuple>

#include <boost/type_index.hpp>

using dbl = double;
using mpfr = long double;


// Description: An interface for all nodes in a function tree, and for a function object as well.  Almost all
//          methods that will be called on a node must be declared in this class.  The main evaluation method is
//          defined in this class.
class Node
{
public:
    ///////// PUBLIC PURE METHODS /////////////////
    
    
    // This method is used to print out a function.
    // If node is an operator, it prints the operation and calls PrintNode on all
    //  children.
    // If node is a symbol, prints out the correct string corresponding to that symbol
    //
    // Every concrete node will implement this method.
    virtual std::string PrintNode() = 0;
    
    // This method adds a child to an operator.
    //
    // Every interface for operator types(nary, binary and unary) implement this method.
    virtual void AddChild(std::shared_ptr<Node> child) = 0;
    
    ///////// PUBLIC PURE METHODS /////////////////
    
    
    
    
    
    // Evaluate the node.  If flag false, just return value, if flag true
    //  run the specific FreshEval of the node, then set flag to false.
    //
    // Template type is type of value you want returned.
    template<typename T>
    T Eval()
    {
        auto& val_pair = std::get< std::pair<T,bool> >(current_value_);
        if(!val_pair.second)
        {
            //            std::cout << "Fresh Eval\n";
            T input{};
            val_pair.first = FreshEval(input);
            val_pair.second = true;
            return val_pair.first;
        }
        else
        {
            return val_pair.first;
        }
    }
    
    
    
    
    
    // TODO(JBC): Should we have a method to Reset all values to fresh Eval?
    
    
    
    // Tells code to run a fresh eval on node for type T
    template <typename T>
    void Reset()
    {
        std::get< std::pair<T,bool> >(current_value_).second = false;
    }

    
    ////////////TESTING////////////////
    static int tabcount;
    virtual void PrintTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << tabcount+1 << "." << boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
    }
    
    ////////////TESTING////////////////

    
    
    
    
protected:
    //Stores the current value of the node in all required types
    //We must hard code in all types that we want here.
    //TODO: Initialize this to some default value, second = false
    std::tuple< std::pair<dbl,bool>, std::pair<mpfr,bool> > current_value_;
    
    
    
    ///////// PRIVATE PURE METHODS /////////////////
    virtual dbl FreshEval(dbl) = 0;
    virtual mpfr FreshEval(mpfr) = 0;
    
    ///////// PRIVATE PURE METHODS /////////////////
    
    
    
    

};



#endif /* defined(__b2Test__Node__) */
