//This file is part of Bertini 2.0.
//
//jacobian.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//jacobian.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with jacobian.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 6/11/15.
//
//
// jacobian.hpp:  Declares the class Jacobian.

/**
\file jacobian.hpp

\brief Provides the Jacobian node type.

*/


#ifndef BERTINI_JACOBIAN_NODE_HPP
#define BERTINI_JACOBIAN_NODE_HPP


#include "bertini2/function_tree/node.hpp"
#include "bertini2/function_tree/symbols/variable.hpp"



namespace bertini {
namespace node{    
		
		/**
		\brief Defines the entry point into a Jacobian tree.

		 This class defines a Jacobian tree.  It stores the entry node for a particular functions tree.
		 */
		class Jacobian : public virtual Function
		{
			friend detail::FreshEvalSelector<dbl>;
			friend detail::FreshEvalSelector<mpfr>;
		public:
				
				
				
				
				
				/**
				 Constructor defines entry node at construct time.
				 */
				Jacobian(const std::shared_ptr<Node> & entry) : Function(entry)
				{
				}
				
				
				/**
				 Jacobians must be evaluated with EvalJ, so that when current_diff_variable changes
				 the Jacobian is reevaluated.
				 */
				template<typename T>
				T Eval(std::shared_ptr<Variable> const& diff_variable = nullptr) const = delete;
				
				
				// Evaluate the node.  If flag false, just return value, if flag true
				//  run the specific FreshEval of the node, then set flag to false.
				template<typename T>
				T EvalJ(std::shared_ptr<Variable> const& diff_variable) const
				{
						auto& val_pair = std::get< std::pair<T,bool> >(current_value_);

						if(diff_variable == current_diff_variable_ && std::get< std::pair<T,bool> >(current_value_).second)
							return val_pair.first;
						else
						{
							current_diff_variable_ = diff_variable;
							Reset();
							std::get< std::pair<T,bool> >(current_value_).first  = detail::FreshEvalSelector<T>::Run(*this,diff_variable);
							std::get< std::pair<T,bool> >(current_value_).second = true;
							return std::get< std::pair<T,bool> >(current_value_).first;
						}						
				}
				



				/**
				 The function which flips the fresh eval bit back to fresh.
				 */
				void Reset() const override
				{
					EnsureNotEmpty();
					
					Node::ResetStoredValues();
					entry_node_->Reset();
				}

				
				virtual ~Jacobian() = default;
				
	
				mutable std::shared_ptr<Variable> current_diff_variable_;
				Jacobian() = default;
		private:
				/**
				 The default constructor
				 */
				
				friend class boost::serialization::access;
		
				template <typename Archive>
				void serialize(Archive& ar, const unsigned version) {
						ar & boost::serialization::base_object<Function>(*this);
				}
		};
		
} // re: namespace node
} // re: namespace bertini



#endif
