//This file is part of Bertini 2.0.
//
//node.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//node.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with node.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// node.hpp:  Declares the class Node.


#ifndef __b2Test__Node__
#define __b2Test__Node__




#include <iostream>
#include <string>
#include <tuple>

#include <boost/type_index.hpp>

#include <complex>
#include "mpfr_complex.hpp"
using dbl = std::complex<double>;
using mpfr = bertini::complex;





namespace bertini {
class Variable;


// Description: An interface for all nodes in a function tree, and for a function object as well.  Almost all
//          methods that will be called on a node must be declared in this class.  The main evaluation method is
//          defined in this class.
class Node
{
public:
	
	virtual ~Node() = default;
	
    ///////// PUBLIC PURE METHODS /////////////////
    
    
    // This method is used to print out a function.
    // If node is an operator, it prints the operation and calls PrintNode on all
    //  children.
    // If node is a symbol, prints out the correct string corresponding to that symbol
    //
    // Every concrete node will implement this method.
    virtual std::string PrintNode() = 0;
    
    
	
	// Tells code to run a fresh eval on node for type T
	
	virtual void Reset()
	{
		ResetStoredValues();
	};
	
    ///////// END PUBLIC PURE METHODS /////////////////
    
    
    
    
    
    // Evaluate the node.  If flag false, just return value, if flag true
    //  run the specific FreshEval of the node, then set flag to false.
    //
    // Template type is type of value you want returned.
    /**
     Evaluate the jacobian with respect to a particular variable.  If flag false, just return value, if flag true
     run the specific FreshEval of the node, then set flag to false.
     
     Template type is type of value you want returned.
     */
    template<typename T>
    T Eval(std::shared_ptr<Variable> diff_variable = nullptr)
    {
        auto& val_pair = std::get< std::pair<T,bool> >(current_value_);
        if(!val_pair.second)
        {
            T input{};
            val_pair.first = FreshEval(input, diff_variable);
            val_pair.second = true;
        }
 
		
		return val_pair.first;
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

	
	
	///////// PUBLIC PURE METHODS /////////////////
	virtual void print(std::ostream& target) const = 0;
    
    virtual std::shared_ptr<Node> Differentiate() = 0;



    /**
	Compute the degree, optionally with respect to a single variable.
    */
    virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const = 0;



    /**
	 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
    */
	std::vector<int> MultiDegree(std::vector< std::shared_ptr<Variable> > const& vars) const
	{
		
		std::vector<int> deg(vars.size());
		for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
		{
			*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
		}
		return deg;
	}


	int Degree(std::vector< std::shared_ptr<Variable > > const& vars) const
	{
		auto multideg = MultiDegree(vars);
		auto deg = 0;
		std::for_each(multideg.begin(),multideg.end(),[&](int n){
                        deg += n;
 						});
		return deg;
	}

    ///////// PUBLIC PURE METHODS /////////////////
    
protected:
    //Stores the current value of the node in all required types
    //We must hard code in all types that we want here.
    //TODO: Initialize this to some default value, second = false
    std::tuple< std::pair<dbl,bool>, std::pair<mpfr,bool> > current_value_;
    
    
    
    ///////// PRIVATE PURE METHODS /////////////////
//    virtual dbl FreshEval(dbl) = 0;
//    virtual mpfr FreshEval(mpfr) = 0;
    
    virtual dbl FreshEval(dbl, std::shared_ptr<Variable>) = 0;
    virtual mpfr FreshEval(mpfr, std::shared_ptr<Variable>) = 0;
	
	
	
    ///////// END PRIVATE PURE METHODS /////////////////
    
    
    
	void ResetStoredValues()
	{
		std::get< std::pair<dbl,bool> >(current_value_).second = false;
		std::get< std::pair<mpfr,bool> >(current_value_).second = false;
	}

};

	
} // re: namespace bertini



namespace  {
	/**
	 a single layer of indirection, though which to call the overridden virtual print() method which must be defined for each non-abstract Node type.
	 */
	std::ostream& operator<<(std::ostream & out, const bertini::Node& N)
	{
		N.print(out);
		return out;
	}
} // re: namespace {}



#endif /* defined(__b2Test__Node__) */






