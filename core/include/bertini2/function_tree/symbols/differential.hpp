//This file is part of Bertini 2.0.
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
#include "function_tree/node.hpp"
#include "function_tree/symbols/symbol.hpp"
#include "function_tree/symbols/number.hpp"



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
		Differential(std::shared_ptr<Variable> diff_variable, std::string var_name) : NamedSymbol(var_name)
		{
			differential_variable_ = diff_variable;
			name("d" + var_name);
		}





		std::shared_ptr<Variable> GetVariable() const 
		{
			return differential_variable_;
		}


		void print(std::ostream & target) const override
		{
			target << name();
		}



		/**
		 Differentiates a number.  Should this return the special number Zero?
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return std::make_shared<Integer>(0);
		}



		virtual ~Differential() = default;




		/**
		Compute the degree with respect to a single variable.   For differentials, the degree is 0.
		*/
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return 0;
		}

		int Degree(VariableGroup const& vars) const override
		{
			return 0;
		}
		
		std::vector<int> MultiDegree(VariableGroup const& vars) const override
		{
			return std::vector<int>(vars.size(),0);
		}

		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			
		}
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return true;
		}

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override
		{
			return true;
		}
		

		/**
		 Change the precision of this variable-precision tree node.
		 
		 \param prec the number of digits to change precision to.
		 */
		virtual void precision(unsigned int prec) override
		{
			auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
			val_pair.first.precision(prec);
		}

		
	protected:
		// This should never be called for a Differential.  Only for Jacobians.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			if(differential_variable_ == diff_variable)
			{
				return 1.0;
			}
			else
			{
				return 0.0;
			}
		}

		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			if(differential_variable_ == diff_variable)
			{
				return mpfr(1);
			}
			else
			{
				return mpfr(0);
			}
		}


	private:
		Differential() = default;
		std::shared_ptr<Variable> differential_variable_;

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<NamedSymbol>(*this);
			ar & differential_variable_;
		}

	};

} // re: namespace node
} // re: namespace bertini

#endif
