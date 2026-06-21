//This file is part of Bertini 2.
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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
//
// silviana amethyst
// University of Wisconsin - Eau Claire
//
//  Created by Collins, James B. on 4/30/15.


/**
\file node.hpp

\brief Defines the abstract Node base class.

*/

#ifndef BERTINI_NODE_BASE_HPP
#define BERTINI_NODE_BASE_HPP




#include <iostream>
#include <string>
#include <tuple>
#include <optional>
#include <typeinfo>

#include <boost/type_index.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/detail/visitable.hpp"
#include "bertini2/function_tree/forward_declares.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>





namespace bertini {

	namespace node{
		class Variable;
	}

using VariableGroup = std::vector< std::shared_ptr<node::Variable> >;



enum class VariableGroupType
{
	Homogeneous,
	Affine,
	Ungrouped
};


namespace node{

namespace detail{
	template<typename T>
	struct FreshEvalSelector
	{};

	template<>
	struct FreshEvalSelector<dbl>
	{
		
		template<typename N>
		static dbl Run(N const& n, std::shared_ptr<Variable> const& diff_variable)
		{
			return n.FreshEval_d(diff_variable);
		}
		
		
		template<typename N>
		static void RunInPlace(dbl& evaluation_value, N const& n, std::shared_ptr<Variable> const& diff_variable)
		{
			n.FreshEval_d(evaluation_value, diff_variable);
		}

	};

	template<>
	struct FreshEvalSelector<mpfr_complex>
	{
		template<typename N>
		static mpfr_complex Run(N const& n, std::shared_ptr<Variable> const& diff_variable)
		{
			return n.FreshEval_mp(diff_variable);
		}
		
		
		template<typename N>
		static void RunInPlace(mpfr_complex& evaluation_value, N const& n, std::shared_ptr<Variable> const& diff_variable)
		{
			n.FreshEval_mp(evaluation_value, diff_variable);
		}

	};
}

/**
\brief Operator precedence classes, used to decide parenthesization when printing.

Higher binds tighter.  A child is wrapped in parentheses only when its
precedence is too low for the position it is printed in; leaves and
self-delimiting nodes (function calls like sin(...), complex pairs) are
atoms and never wrapped.  Negative literal constants report PrecNegate so
they parenthesize exactly where a Negate node would.
*/
enum PrintPrecedence : unsigned {
	PrecSum = 10,
	PrecNegate = 15,
	PrecMult = 20,
	PrecPower = 30,
	PrecAtom = 100
};

/**
An interface for all nodes in a function tree, and for a function object as well.  Almost all
 methods that will be called on a node must be declared in this class.  The main evaluation method is
 defined in this class.

 Bertini function trees are created using std::shared_ptr's to Node base class, generally. 

 \brief Abstract base class for the Bertini hybrid-precision (double-multiple) expression tree. 
 */
class Node : public VisitableBase<>, public std::enable_shared_from_this<Node>
{
	friend detail::FreshEvalSelector<dbl>;
	friend detail::FreshEvalSelector<mpfr_complex>;
public:
	
	virtual ~Node() = default;
	
	///////// PUBLIC PURE METHODS /////////////////
  


	/**
	 Tells code to run a fresh eval on node next time.
	*/
	virtual void Reset() const = 0;
	
	///////// END PUBLIC PURE METHODS /////////////////
	
	

	public:
	/**
	 Evaluate the node.  If flag false, just return value, if flag true
	 run the specific FreshEval of the node, then set flag to false.

	 Template type is type of value you want returned.

	 \return The value of the node.
	 \tparam T The number type for return.  Must be one of the types stored in the Node class, currently dbl and mpfr_complex.
	 */
	template<typename T>
	T Eval(std::shared_ptr<Variable> const& diff_variable = nullptr) const 
	{
		T result;
		EvalInPlace(result, diff_variable);
		return result;
	}
	

	/**
	 Evaluate the node in place.  If flag false, just return value, if flag true
	 run the specific FreshEval of the node, then set flag to false.
	 
	 Template type is type of value you want returned.
	 
	 \return The value of the node.
	 \tparam T The number type for return.  Must be one of the types stored in the Node class, currently dbl and mpfr_complex.
	 */
	template<typename T>
	void EvalInPlace(T& eval_value, std::shared_ptr<Variable> const& diff_variable = nullptr) const;

	
	///////// PUBLIC PURE METHODS /////////////////



	/**
	\brief Functionally simplify this tree, returning a NEW simplified tree.

	Non-mutating successor to the in-place EliminateZeros / EliminateOnes / ReduceDepth
	machinery (ADR-0011, issue 251).  The input tree is never modified; a fresh simplified
	tree is returned, built through the SimplifiedSum / SimplifiedMult / SimplifiedNegate
	factories so that literal zeros/ones vanish and exact constants fold.  The default
	(leaves, and any node with nothing to simplify) returns the node unchanged -- structural
	sharing preserved.

	\return A simplified tree (possibly the same node, if nothing simplified).
	*/
	virtual std::shared_ptr<Node> Simplified() const;

	/**
	Virtual method for printing Nodes to arbitrary output streams.
	*/
	virtual void print(std::ostream& target) const = 0;

	/**
	\brief The printing precedence of this node, deciding parenthesization.

	Defaults to PrecAtom: leaves and self-delimiting nodes are never wrapped.
	Operator nodes override this; printing parents wrap a child only when its
	precedence is too low for the position it occupies.
	*/
	virtual unsigned Precedence() const
	{
		return PrecAtom;
	}

	/**
	\brief Is this node the literal constant 0?

	Exact check on stored literal values (Integer/Float/Rational override);
	false for everything else, including non-literal expressions that happen
	to evaluate to zero.  Used by differentiation to build already-simplified
	trees without evaluating anything.
	*/
	virtual bool IsLiteralZero() const
	{
		return false;
	}

	/**
	\brief Is this node the literal constant 1?

	\see IsLiteralZero
	*/
	virtual bool IsLiteralOne() const
	{
		return false;
	}

	/**
	\brief Memoized structural hash of this node.

	Order-sensitive: equal structures hash equal, where children are folded in by their own
	Hash() and operand order / signs / exponents participate.  The default (leaves and any
	node not overriding HashImpl) is node identity (the object's address), so distinct objects
	hash distinctly; value nodes (Integer/Rational/Float) and operators override to be
	structural.  Memoized -- valid because nodes are immutable post-construction -- and
	deliberately independent of the mutable working precision (stable across precision()).

	This is the predicate layer for the hash-consing intern table; nothing wires it into
	construction yet.
	*/
	std::size_t Hash() const;

	/**
	\brief Order-sensitive structural equality, shallow given interned children.

	Two nodes are the same iff they have the same dynamic type, the same operator payload
	(signs / mult-or-div flags / integer exponent / literal value), and the same operands
	**by pointer** (operands are not recursed -- in the hash-consed world children are already
	canonical, so pointer-equality is structural equality).  The default is node identity;
	value/operator nodes override.  Consistent with Hash(): IsSame(a,b) implies a.Hash()==b.Hash().
	*/
	virtual bool IsSame(Node const& other) const;


	/**
	\brief Compute the derivative with respect to a single variable.

	Virtual method for differentiating the node.  If no variable is passed, produces a Jacobian tree when all is said and done, which is used for evaluating the Jacobian.

	\return The result of differentiation.  Jacobian or regular Node depending on what you passed in.
	*/
	virtual std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const = 0;
	
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
	Homogenize a tree, returning a NEW homogenized tree (functional / non-mutating).  Input a
	variable group holding the non-homogeneous variables, and the new homogenizing variable.
	The homvar may be an element of the variable group, that's perfectly ok.

	The input tree is never modified; degree-deficient summands are padded with powers of homvar
	in a freshly-built tree (so a throw on a non-polynomial term can't leave a half-homogenized
	tree behind).  The default (leaves, and anything with nothing to homogenize) returns the node
	unchanged -- shared variables are preserved, so set_current_value still drives the result.

	\param homvar The homogenizing variable, which is multiplied against terms with degree deficiency with repect to other terms.
	\param vars A group of variables, with respect to which you wish to homogenize.
	\return A homogenized tree (possibly the same node, if nothing changed).
	*/
	virtual std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const;

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
	virtual void precision(unsigned int prec) const = 0;

	unsigned precision() const;

	///////// PUBLIC PURE METHODS /////////////////

	/**
	Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

	\return True if it is polynomial, false if not.
	*/
	bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const;
	

	/**
	Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

	\return True if it is polynomial, false if not.
	*/
	bool IsPolynomial(VariableGroup const&v) const;


	
	
	

protected:
	//Stores the current value of the node in all required types
	//We must hard code in all types that we want here.
	//TODO: Initialize this to some default value, second = false
	mutable std::tuple< std::pair<dbl,bool>, std::pair<mpfr_complex,bool> > current_value_;

	/// Memoized structural hash (computed on first Hash() call; nodes are immutable so it
	/// never needs invalidating).  Not serialized -- a clone recomputes it on demand.
	mutable std::optional<std::size_t> structural_hash_;

	/// Compute this node's structural hash.  Default: identity (the object address), so
	/// distinct nodes hash distinctly.  Value/operator nodes override to be structural.
	virtual std::size_t HashImpl() const;

	/// Combine an extra value into a running hash (boost::hash_combine recipe).
	static void HashCombine(std::size_t& seed, std::size_t value)
	{
		seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
	}



	///////// PRIVATE PURE METHODS /////////////////
	
	/**
	Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper Eval<>() call from Node, if so required (by resetting, etc).

	If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr_complex.
	*/
	virtual dbl FreshEval_d(std::shared_ptr<Variable> const&) const = 0;

	/**
	 Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper EvalInPlace<>() call from Node, if so required (by resetting, etc).
	 
	 If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr_complex.
	 */
	virtual void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const&) const = 0;

	
	/**
	Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper Eval<>() call from Node, if so required (by resetting, etc).

	If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr_complex.
	*/
	virtual mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const&) const = 0;
	
	/**
	 Overridden code for specific node types, for how to evaluate themselves.  Called from the wrapper Eval<>() call from Node, if so required (by resetting, etc).
	 
	 If we had the ability to use template virtual functions, we would have.  However, this is impossible with current C++ without using experimental libraries, so we have two copies -- because there are two number types for Nodes, dbl and mpfr_complex.
	 */
	virtual void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const&) const = 0;

	
	
	
	///////// END PRIVATE PURE METHODS /////////////////
	
	
	/**
	Set the stored values for the Node to indicate a fresh eval on the next pass.  This is so that Nodes which are referred to more than once, are only evaluated once.  The first evaluation is fresh, and then the indicator for fresh/stored is set to stored.  Subsequent evaluation calls simply return the stored number.
	*/
	void ResetStoredValues() const;

	Node();

private:
	friend std::ostream& operator<<(std::ostream & out, const Node& N);

	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned /*version*/) {
		register_derived_node_types(ar);
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
	
	inline std::ostream& operator<<(std::ostream & out, const std::shared_ptr<Node>& N)
	{
		N->print(out);
		return out;
	}



	} // re: namespace node
} // re: namespace bertini



#endif 
/* defined(BERTINI_NODE_BASE_HPP) */






