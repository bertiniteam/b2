//This file is part of Bertini 2.
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
#include "bertini2/function_tree/roots/function.hpp"
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
				 */
				Jacobian(const std::shared_ptr<Node> & entry);
				
				
				/**
				 Jacobians must be evaluated with EvalJ, so that when current_diff_variable changes
				 the Jacobian is reevaluated.
				 */
				template<typename T>
				T Eval(std::shared_ptr<Variable> const& diff_variable = nullptr) const = delete;
				
				
				// Evaluate the node.  If flag false, just return value, if flag true
				//  run the specific FreshEval of the node, then set flag to false.
				template<typename T>
				T EvalJ(std::shared_ptr<Variable> const& diff_variable) const;				


				// Evaluate the node.  If flag false, just return value, if flag true
				//  run the specific FreshEval of the node, then set flag to false.
				template<typename T>
				void EvalJInPlace(T& eval_value, std::shared_ptr<Variable> const& diff_variable) const;


				/**
				 The function which flips the fresh eval bit back to fresh.
				 */
				void Reset() const override;

				
				virtual ~Jacobian() = default;
				
				/**
				\brief Default construction of a Jacobian node is forbidden
				*/
				Jacobian() = default;

		private:

			mutable std::shared_ptr<Variable> current_diff_variable_;


			friend class boost::serialization::access;
	
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
					ar & boost::serialization::base_object<Function>(*this);
			}
		};
		
} // re: namespace node
} // re: namespace bertini



#endif
