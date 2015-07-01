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


#ifndef b2Test_Function_h
#define b2Test_Function_h


#include "function_tree/symbols/symbol.hpp"



namespace bertini {
	
	
	/**
	 Node -> Function
	 This class defines a function.  It stores the entry node for a particular functions tree.
	 */
	class Function : public NamedSymbol
	{
	public:
		
		
		/**
		 The default constructor
		 */
		Function() {};
		
		

		
		Function(std::string new_name) : NamedSymbol(new_name)
		{ }
		
		
		/**
		 Constructor defines entry node at construct time.
		 */
		Function(const std::shared_ptr<Node> & entry) : entry_node_(entry)
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
		 Virtual overridden function for piping the tree to an output stream.
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
        virtual std::shared_ptr<Node> Differentiate() override
        {
            return entry_node_->Differentiate();
        }
		
		


		/**
        Compute the degree of a node.  For functions, the degree is the degree of the entry node.
        */
        virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
        {
            return entry_node_->Degree(v);
        }


        virtual void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			entry_node_->Homogenize(vars, homvar);
		}


		
	protected:
		
		/**
         Calls FreshEval on the entry node to the tree.
         */
        virtual dbl FreshEval(dbl d, std::shared_ptr<Variable> diff_variable) override
        {
            return entry_node_->Eval<dbl>(diff_variable);
        }
        
        /**
         Calls FreshEval on the entry node to the tree.
         */
        virtual mpfr FreshEval(mpfr m, std::shared_ptr<Variable> diff_variable) override
        {
            return entry_node_->Eval<mpfr>(diff_variable);
        }

		
		
		std::shared_ptr<Node> entry_node_; ///< The top node for the function.
	};
	
	
} // re: namespace bertini



#endif
