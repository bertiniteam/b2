//This file is part of Bertini 2.
//
//include/bertini2/function_tree/symbols/variable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/symbols/variable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/symbols/variable.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
// include/bertini2/function_tree/symbols/variable.hpp:  Declares the class Variable.


/**
\file include/bertini2/function_tree/symbols/variable.hpp

\brief Provides the Variable Node class.

*/
#ifndef BERTINI_FUNCTION_TREE_VARIABLE_HPP
#define BERTINI_FUNCTION_TREE_VARIABLE_HPP

#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/symbols/differential.hpp"
#include "bertini2/function_tree/factory.hpp"


namespace  bertini {
namespace node{
	/**
	\brief Represents variable leaves in the function tree.

	This class represents variable leaves in the function tree.  FreshEval returns
	the current value of the variable.

	When differentiated, produces a differential referring to it.
	*/
	class Variable : public virtual NamedSymbol, public std::enable_shared_from_this<Variable>
	{
	public:
		
		
		Variable(std::string new_name);
		
		
		
		
		virtual ~Variable() = default;
		


		explicit operator std::string(){return name();}
		
		
		
		// This sets the value for the variable
		template <typename T>
		void set_current_value(T const& val);

		/**
		\brief Changes the value of the variable to be not-a-number.  
		*/
		template <typename T>
		void SetToNan();

		/**
		\brief Changes the value of the variable to be a random complex number.  
		*/
		template <typename T>
		void SetToRand();

		/**
		\brief Changes the value of the variable to be a random complex number, of magnitude 1.  
		*/
		template <typename T>
		void SetToRandUnit();

		/**
		 Differentiates a variable.  
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		void Reset() const override;

		/**
		Compute the degree with respect to a single variable.

		If this is the variable, then the degree is 1.  Otherwise, 0.
		*/
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;


		int Degree(VariableGroup const& vars) const override;

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
		
		// Return current value of the variable.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		Variable();
	private:
		
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<NamedSymbol>(*this);
		}

	};
	

	

} // re: namespace node
} // re: namespace bertini




#endif
