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

using VariableGroup = std::vector< std::shared_ptr<Variable> >;

// Description: An interface for all nodes in a function tree, and for a function object as well.  Almost all
//          methods that will be called on a node must be declared in this class.  The main evaluation method is
//          defined in this class.
class Node
{
public:
	
	virtual ~Node() = default;
	
    ///////// PUBLIC PURE METHODS /////////////////
  
	
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
	virtual std::vector<int> MultiDegree(std::vector< std::shared_ptr<Variable> > const& vars) const = 0;

	/**
	 Compute the overall degree with respect to a variable group.
	*/
	virtual int Degree(std::vector< std::shared_ptr<Variable > > const& vars) const = 0;

	/**
	Homogenize a tree, inputting a variable group holding the non-homogeneous variables, and the new homogenizing variable.  The homvar may be an element of the variable group, that's perfectly ok.
	*/
	virtual void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) = 0;

	/**
	Check for homogeneity, absolutely with respect to all variables, including path variables and all other variable types, or with respect to a single varaible, if passed. 
	*/
	virtual bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const = 0;

	/**
	Check for homogeneity, with respect to a variable group.
	*/
	virtual bool IsHomogeneous(VariableGroup const& vars) const = 0;
    ///////// PUBLIC PURE METHODS /////////////////
    

	bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const
	{
		return Degree(v)>=0;
	}

	bool IsPolynomial(VariableGroup const&v) const
	{
		return Degree(v)>=0;
	}

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
	
	/**
	 a single layer of indirection, though which to call the overridden virtual print() method which must be defined for each non-abstract Node type.
	 */
	inline std::ostream& operator<<(std::ostream & out, const bertini::Node& N)
	{
		N.print(out);
		return out;
	}
	
} // re: namespace bertini



namespace  {
	
} // re: namespace {}



#endif /* defined(__b2Test__Node__) */






