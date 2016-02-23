//This file is part of Bertini 2.0.
//
//variable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//variable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with variable.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// variable.hpp:  Declares the class Variable.


/**
\file variable.hpp

\brief Provides the Variable Node class.

*/
#ifndef BERTINI_FUNCTION_TREE_VARIABLE_HPP
#define BERTINI_FUNCTION_TREE_VARIABLE_HPP

#include "function_tree/symbols/symbol.hpp"
#include "function_tree/symbols/differential.hpp"



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
		
		
		Variable(std::string new_name) : NamedSymbol(new_name)
		{ }
		
		
		
		
		virtual ~Variable() = default;
		


		explicit operator std::string(){return name();}
		
		
		
		// This sets the value for the variable
		template <typename T>
		void set_current_value(T val)
		{
			std::get< std::pair<T,bool> >(current_value_).first = val;
			std::get< std::pair<T,bool> >(current_value_).second = false;
		}
		
		
		/**
		 Differentiates a variable.  
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return std::make_shared<Differential>(shared_from_this(), name());
		}
		
		

		


		/**
		Compute the degree with respect to a single variable.

		If this is the variable, then the degree is 1.  Otherwise, 0.
		*/
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			if (v)
			{
				if (this == v.get())
					return 1;
				else
					return 0;
			}
			else
				return 1;
			
		}


		int Degree(VariableGroup const& vars) const override
		{
			for (auto iter : vars)
				if (this==iter.get())
					return 1;
				
			return 0;
		}

		std::vector<int> MultiDegree(VariableGroup const& vars) const override
		{
			std::vector<int> deg;
			for (auto iter=vars.begin(); iter!=vars.end(); iter++)
				if (this==(*iter).get())
					deg.push_back(1);
				else
					deg.push_back(0);
			return deg;
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
		
		// Return current value of the variable.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return std::get< std::pair<dbl,bool> >(current_value_).first;
		}
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return std::get< std::pair<mpfr,bool> >(current_value_).first;
		}
		Variable() = default;
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
