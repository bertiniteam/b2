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


#include "function_tree/node.hpp"



namespace bertini {
	
	
	/**
	 Node -> Function
	 This class defines a function.  It stores the entry node for a particular functions tree.
	 */
	class Function : public Node
	{
	public:
		
		
		/**
		 The default constructor
		 */
		Function()
		{};
		
		
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
			EnsureNotEmpty();
			entry_node_->print(target);
		}
		
		
		virtual ~Function() = default;
		
		
		
		/**
		 This returns a string that holds a print out of the function.
		 */
		virtual std::string PrintNode() override
		{
			EnsureNotEmpty();
			return entry_node_->PrintNode();
		}
		
		
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
		
		
		////////////// TESTING /////////////////
		virtual void PrintTree() override
		{
			entry_node_->PrintTree();
		}
		////////////// TESTING /////////////////
		
	private:
		
		/**
		 Calls FreshEval on the entry node to the tree.
		 */
		virtual dbl FreshEval(dbl d)
		{
			return entry_node_->Eval<dbl>();
		}
		
		/**
		 Calls FreshEval on the entry node to the tree.
		 */
		virtual mpfr FreshEval(mpfr m)
		{
			return entry_node_->Eval<mpfr>();
		}
		
		
		
		std::shared_ptr<Node> entry_node_; ///< The top node for the function.
	};
	
	
} // re: namespace bertini



#endif
