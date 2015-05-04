//
//  main.cpp
//  Variadic Templates
//
//  Created by Collins, James B. on 4/23/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <iostream>
#include <vector>


/* Abstract Node class.
 Base for all other Nodes */

template<typename... TArgs>
class Node
{
    
public:
    
    /* This is the eval function that actually get's called
     It then calls a particular eval based on what DerivedType is*/
    
    template<typename T, typename DerivedType>
    T eval()
    {
        return static_cast<DerivedType*>(this)->template eval<T>();
    }
    
    
    template<typename T>
    void set(T val)
    {
        std::get< std::pair<T,bool> >(this->current_value).first = val;
    }
    
    
    /* Add a child to the tree */
    
    void addChild(Node<TArgs...>* child)
    {
        children.push_back(child);
    }
    
    virtual auto typeofNode() -> decltype(this);
//    virtual double eval(double) = 0;
//    virtual int eval(int) = 0;
    
    
    
    
protected:
    std::tuple< std::pair<TArgs,bool>... > current_value;
    std::vector< Node<TArgs...>* > children;
};



template<typename... TArgs>
class Value : public Node<TArgs...>
{
public:
    //    Sum()
    //    {
    //        this->current_value = std::pair<float,bool> {7.9,true};
    //    }
    
    
    
    
    template<typename T>
    T eval()
    {
        return std::get< std::pair<T,bool> >(this->current_value).first;
    }
    
    
    
    virtual auto typeofNode() -> decltype(this);
};


template<typename... TArgs>
class Sum : public Node<TArgs...>
{
public:
//    Sum()
//    {
//        this->current_value = std::pair<float,bool> {7.9,true};
//    }
    
    
    
    
    template<typename T>
    T eval()
    {
        T retval = 0;

        for(auto c : this->children)
        {
            auto test = static_cast<Value<TArgs...> *>(c)->typeofNode();
            T child = c->template eval<T, Value<TArgs...> >();  //I don't understand why template is needed, compiler made me do it
            retval += child;
        }
        
        return retval;
    }
};






//template<typename... TArgs>
//class Mult : public Node<Mult<TArgs...>, TArgs...>
//{
//public:
////    Mult()
////    {
////        this->current_value = std::pair<float,bool> {7.9,true};
////    }
//    
//    
//    template<typename T>
//    T eval()
//    {
//        return std::get< std::pair<T,bool> >(this->current_value).first * 2;
//    }
//};




int main(int argc, const char * argv[]) {
    // insert code here...
    
    Node<float, double>* x = new Sum<float, double>();
    Node<float, double>* y1 = new Value<float, double>();
    Node<float, double>* y2 = new Value<float, double>();
    
    // Set values in all nodes
    x->set<float>(1.2f);
    x->set<double>(3.1);
    
    y1->set<float>(.1f);
    y1->set<double>(.3);
    
    y2->set<float>(1.1f);
    y2->set<double>(1.3);
    
    
    //Add the children y1 & y2 to parent x
    x->addChild(y1);
    x->addChild(y2);
    
    std::cout << x->eval<float, Sum<float, double> >() << std::endl;
    std::cout << x->eval<double, Sum<float, double> >() << std::endl;
    
    
    return 0;
}
