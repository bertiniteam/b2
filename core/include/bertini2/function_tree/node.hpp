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

/**
\file node.hpp

\brief Defines the abstract Node base class.

*/

#ifndef BERTINI_NODE_BASE_HPP
#define BERTINI_NODE_BASE_HPP




#include <iostream>
#include <string>
#include <tuple>

#include <boost/type_index.hpp>

#include "num_traits.hpp"


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include <deque>




namespace bertini {

	namespace node{
		class Variable;
	}

using VariableGroup = std::deque< std::shared_ptr<node::Variable> >;

enum class VariableGroupType
{
	Homogeneous,
	Affine,
	Ungrouped
};


namespace node{


/**
An interface for all nodes in a function tree, and for a function object as well.  Almost all
 methods that will be called on a node must be declared in this class.  The main evaluation method is
 defined in this class.

 Bertini function trees are created using std::shared_ptr's to Node base class, generally. 

 \brief Abstract base class for the Bertini hybrid-precision (double-multiple) expression tree. 
 */
class Node
{
public:
	
	virtual ~Node() = default;
	
	///////// PUBLIC PURE METHODS /////////////////
  


	/**
	 Tells code to run a fresh eval on node next time.
	*/
	virtual void Reset()
	{
		ResetStoredValues();
	};
	
	///////// END PUBLIC PURE METHODS /////////////////
	
	
	
   
	/**
	 Evaluate the node.  If flag false, just return value, if flag true
	 run the specific FreshEval of the node, then set flag to false.

	 Template type is type of value you want returned.

	 \return The value of the node.
	 \tparam T The number type for return.  Must be one of the types stored in the Node class, currently dbl and mpfr.
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
	/**
	Virtual method for printing Nodes to arbitrary output streams.
	*/
	virtual void print(std::ostream& target) const = 0;
	

	/**
	Virtual method for differentiating the node.  Produces a Jacobian tree when all is said and done, which is used for evaluating the Jacobian.

	\return The Jacobian for the node.
	*/
	virtual std::shared_ptr<Node> Differentiate() = 0;



	/**
	Compute the degree, optionally with respect to a single variable.

	\param v Shared pointer to variable with respect to which you want to compute the degree of the Node.
	\return The degree.  Will be negative if the Node is non-polynomial.
	*/
	virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const = 0;



	/**
	 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  

	 \return A vector containing the degrees.  Negative entries indicate non-polynomiality.
	*/
	virtual std::vector<int> MultiDegree(VariableGroup const& vars) const = 0;

	/**
	 Compute the overall degree with respect to a variable group.
	
	\param vars A group of variables.
	 \return The degree.  Will be negative if the Node is non-polynomial.
	*/
	virtual int Degree(VariableGroup const& vars) const = 0;

	/**
	Homogenize a tree, inputting a variable group holding the non-homogeneous variables, and the new homogenizing variable.  The homvar may be an element of the variable group, that's perfectly ok.
	
	\param homvar The homogenizing variable, which is multiplied against terms with degree deficiency with repect to other terms.
	\param vars A group of variables, with respect to which you wish to homogenize.
	*/
	virtual void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) = 0;

	/**
	Check for homogeneity, absolutely with respect to all variables, including path variables and all other variable types, or with respect to a single varaible, if passed. 
	
	\return True if it is homogeneous, false if not.
	*/
	virtual bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const = 0;

	/**
	Check for homogeneity, with respect to a variable group.

	\return True if it is homogeneous, false if not.
	*/
	virtual bool IsHomogeneous(VariableGroup const& vars) const = 0;
	
	
	/**
	 Change the precision of this variable-precision tree node.
	 
	 \param prec the number of digits to change precision to.
	 */
	virtual void precision(unsigned int prec) = 0;

	///////// PUBLIC PURE METHODS /////////////////

	/**
	Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

	\return True if it is polynomial, false if not.
	*/
	bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const
	{
		return Degree(v)>=0;
	}


	/**
	Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

	\return True if it is polynomial, false if not.
	*/
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
	
	/**
	Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper Eval<>() call from Node, if so required (by resetting, etc).

	If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr.
	*/
	virtual dbl FreshEval(dbl, std::shared_ptr<Variable>) = 0;

	/**
	Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper Eval<>() call from Node, if so required (by resetting, etc).

	If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr.
	*/
	virtual mpfr FreshEval(mpfr, std::shared_ptr<Variable>) = 0;
	
	
	
	///////// END PRIVATE PURE METHODS /////////////////
	
	
	/**
	Set the stored values for the Node to indicate a fresh eval on the next pass.  This is so that Nodes which are referred to more than once, are only evaluated once.  The first evaluation is fresh, and then the indicator for fresh/stored is set to stored.  Subsequent evaluation calls simply return the stored number.
	*/
	void ResetStoredValues()
	{
		std::get< std::pair<dbl,bool> >(current_value_).second = false;
		std::get< std::pair<mpfr,bool> >(current_value_).second = false;
	}

private:
	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned version) {
	}

}; // re: class node
	
	/**
	 a single layer of indirection, though which to call the overridden virtual print() method which must be defined for each non-abstract Node type.
	 */
	inline std::ostream& operator<<(std::ostream & out, const Node& N)
	{
		N.print(out);
		return out;
	}
	
	} // re: namespace node
} // re: namespace bertini







#endif 
/* defined(BERTINI_NODE_BASE_HPP) */






