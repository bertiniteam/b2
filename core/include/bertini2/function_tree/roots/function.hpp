//This file is part of Bertini 2.0.
//
//function.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// function.hpp:  Declares the class Function.

/**
\file function.hpp

\brief Provides the Function Node type, a NamedSymbol.

*/


#ifndef BERTINI_FUNCTION_NODE_HPP
#define BERTINI_FUNCTION_NODE_HPP


#include "function_tree/symbols/symbol.hpp"



namespace bertini {
namespace node{
	
	/**
	\brief Formal entry point into an expression tree.

	This class defines a function.  It stores the entry node for a particular functions tree.
	 */
	class Function : public NamedSymbol
	{
	public:
		
		
		
		

		
		Function(std::string new_name) : NamedSymbol(new_name)
		{ }
		
		
		/**
		 Constructor defines entry node at construct time.
		 */
		Function(const std::shared_ptr<Node> & entry) : NamedSymbol("unnamed_function"), entry_node_(entry)
		{
		}
		
		
		
		/**
		 Get the pointer to the entry node for this function.
		 */
		std::shared_ptr<Node> entry_node() const
		{
			return entry_node_;
		}
		
		
		/**
		 overridden function for piping the tree to an output stream.
		 */
		void print(std::ostream & target) const override
		{
			if (entry_node_)
				entry_node_->print(target);
			else
				target << "emptyfunction";
		}
		
		
		virtual ~Function() = default;
		
		
		
		
		
		/**
		 The function which flips the fresh eval bit back to fresh.
		 */
		void Reset() override
		{
			EnsureNotEmpty();
			
			Node::ResetStoredValues();
			entry_node_->Reset();
		}
		
		
		/**
		 Add a child onto the container for this operator
		 */
		void SetRoot(std::shared_ptr<Node> entry)
		{
			entry_node_ = entry;
		}
		
		
		/**
		 throws a runtime error if the root node is nullptr
		 */
		void EnsureNotEmpty() const
		{
			if (entry_node_==nullptr)
			{
				throw std::runtime_error("Function node type has empty root node");
			}
		}
		
		
		/** 
		 Calls Differentiate on the entry node and returns differentiated entry node.
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return entry_node_->Differentiate();
		}
		
		


		/**
		Compute the degree of a node.  For functions, the degree is the degree of the entry node.
		*/
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return entry_node_->Degree(v);
		}


		int Degree(VariableGroup const& vars) const override
		{
			return entry_node_->Degree(vars);
		}



		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
		*/
		std::vector<int> MultiDegree(VariableGroup const& vars) const override
		{
			
			std::vector<int> deg(vars.size());
			for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
			{
				*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
			}
			return deg;
		}


		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			entry_node_->Homogenize(vars, homvar);
		}

		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return entry_node_->IsHomogeneous(v);
		}

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override
		{
			return entry_node_->IsHomogeneous(vars);
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
		
		/**
		 Calls FreshEval on the entry node to the tree.
		 */
		dbl FreshEval(dbl d, std::shared_ptr<Variable> diff_variable) override
		{
			return entry_node_->Eval<dbl>(diff_variable);
		}
		
		/**
		 Calls FreshEval on the entry node to the tree.
		 */
		mpfr FreshEval(mpfr m, std::shared_ptr<Variable> diff_variable) override
		{
			return entry_node_->Eval<mpfr>(diff_variable);
		}

		
		
		std::shared_ptr<Node> entry_node_; ///< The top node for the function.
		
		Function() = default;
	private:
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<NamedSymbol>(*this);
			ar & entry_node_;
		}
	};
	
} // re: namespace node
} // re: namespace bertini



#endif
