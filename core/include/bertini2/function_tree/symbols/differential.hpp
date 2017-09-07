//This file is part of Bertini 2.
//
//differential.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//differential.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with differential.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
// differential.hpp:  Declares the class SpecialNumber.

/**
\file differential.hpp

\brief Provides the differential Node class.

*/

#ifndef BERTINI_NODE_DIFFERENTIAL_HPP
#define BERTINI_NODE_DIFFERENTIAL_HPP

#include <memory>
#include "bertini2/function_tree/node.hpp"
#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/symbols/number.hpp"



namespace bertini {
namespace node{
	class Variable;


	/**
	\brief Provides the differential type for differentiation of expression trees.

	This class represents differentials.  These are produced in the course of differentiation of a non-constant expression tree.
	*/
	class Differential : public virtual NamedSymbol
	{
	public:
		


		/**
		 Input shared_ptr to a Variable.
		 */
		Differential(std::shared_ptr<const Variable> diff_variable, std::string var_name);


		void Reset() const override;


		const std::shared_ptr<const Variable>& GetVariable() const;


		void print(std::ostream & target) const override;



		/**
		 Differentiates a number.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~Differential() = default;




		/**
		Compute the degree with respect to a single variable.   For differentials, the degree is 0.
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
		// This should never be called for a Differential.  Only for Jacobians.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;



	private:
		Differential() = default;
		std::shared_ptr<const Variable> differential_variable_;


		friend class boost::serialization::access;
		
		template<class Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			ar & boost::serialization::base_object<NamedSymbol>(*this);
			ar & differential_variable_;
			// ar & const_cast<std::shared_ptr<const Variable> >(differential_variable_);
		}
		
		template<class Archive>
		void load(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<NamedSymbol>(*this);
			ar & std::const_pointer_cast<Variable>(differential_variable_);
		}
		
		// BOOST_SERIALIZATION_SPLIT_MEMBER()

	};

} // re: namespace node
} // re: namespace bertini

#endif
