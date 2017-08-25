//This file is part of Bertini 2.
//
//operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// Dani Brake
// University of Notre Dame
//
//  Created by Collins, James B. on 4/30/15.
//
//
// operator.hpp:  Declares the class Operator.

/**
\file operator.hpp

\brief Provides the abstract Operator Node type.

*/

#ifndef BERTINI_OPERATOR_HPP
#define BERTINI_OPERATOR_HPP

#include "bertini2/function_tree/node.hpp"

#include <vector>


namespace bertini {
namespace node{
	
	/**
	\brief Abstract Node type from which all Operators inherit.

	This class is an interface for all operators in the Bertini2 function tree.
	*/
	class Operator : public virtual Node
	{
		
	public:
		
		virtual ~Operator() = default;

	private:

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Node>(*this);
		}
	};
	
	
	
	/**
	\brief Interface for all unary Operator types.


	This class is an interface for all unary operators, such as negation.
	The sole child is stored in a shared_ptr.
	*/
	class UnaryOperator : public virtual Operator
	{
	public:
		
		UnaryOperator(const std::shared_ptr<Node> & N) : child_(N)
		{}
		
		
		
		virtual ~UnaryOperator() = default;
		
		
		void Reset() const override;
		
		
		void SetChild(std::shared_ptr<Node> new_child);
		
		
		
		//Return the only child for the unary operator
		std::shared_ptr<Node> first_child() const;
		
		
		
		/**
		 Compute the degree of a node.  For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return nan.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		int Degree(VariableGroup const& vars) const override;


		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
		*/
		std::vector<int> MultiDegree(VariableGroup const& vars) const override;


		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override;

		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		
		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override;
		


		/**
		 Change the precision of this variable-precision tree node.
		 
		 \param prec the number of digits to change precision to.
		 */
		void precision(unsigned int prec) const override;



	protected:
		//Stores the single child of the unary operator
		std::shared_ptr<Node> child_;
		UnaryOperator(){}
	private:
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Operator>(*this);
			ar & child_;
		}
	};
	
	

	
	/**
	\brief Abstract interface for n-ary Operator types.

	This class is an interface for all n-ary operators, such as summation and multiplication.
	Children of the operator are stored in a vector and methods to add and access children are available
	in this interface.
	*/
	class NaryOperator : public virtual Operator
	{
	public:
		
		virtual ~NaryOperator() = default;
		
		
		void Reset() const override;
		
		// Add a child onto the container for this operator
		virtual void AddChild(std::shared_ptr<Node> child);
		
		
		size_t children_size() const;
		
		std::shared_ptr<Node> first_child() const;
		
		
		
		
		 /**
		 Change the precision of this variable-precision tree node.
		 
		 \param prec the number of digits to change precision to.
		 */
		void precision(unsigned int prec) const override;

		
		
	protected:
		
		//Stores all children for this operator node.
		//This is an NaryOperator and can have any number of children.
		std::vector< std::shared_ptr<Node> > children_;

		// constructor is protected to help prevent instantiating empty operators
		NaryOperator(){}
		
	private:

		virtual void PrecisionChangeSpecific(unsigned prec) const;


		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Operator>(*this);
			ar & children_;
		}
		
	};
	
	
} // re: namespace node	
} // re: namespace bertini


#endif
