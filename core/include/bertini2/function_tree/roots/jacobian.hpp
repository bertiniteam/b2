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


#include "function_tree/node.hpp"
#include "function_tree/symbols/variable.hpp"



namespace bertini {
namespace node{    
		
		/**
		\brief Defines the entry point into a Jacobian tree.

		 This class defines a Jacobian tree.  It stores the entry node for a particular functions tree.
		 */
		class Jacobian : public Function
		{
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
				T Eval(std::shared_ptr<Variable> diff_variable = nullptr) = delete;
				
				
				// Evaluate the node.  If flag false, just return value, if flag true
				//  run the specific FreshEval of the node, then set flag to false.
				//
				// Template type is type of value you want returned.
				template<typename T>
				T EvalJ(std::shared_ptr<Variable> diff_variable)
				{
						//		this->print(std::cout);
						//		std::cout << " has value ";
						
						auto& val_pair = std::get< std::pair<T,bool> >(current_value_);
						auto& variable_pair = std::get< std::pair<T,std::shared_ptr<Variable>> >(current_diff_variable_);
						if(diff_variable != variable_pair.second)
						{
								//            std::cout << "Fresh Eval\n";
								variable_pair.second = diff_variable;
								Reset();
								T input{};
								val_pair.first = FreshEval(input, diff_variable);
								val_pair.second = true;
						}
						
						
						//		std::cout << val_pair.first << std::endl;
						return val_pair.first;
				}
				



				

				
				virtual ~Jacobian() = default;
				
			
				
				
	//       /**
		//  Change the precision of this variable-precision tree node.
		 
		//  \param prec the number of digits to change precision to.
		//  */
		// virtual void precision(unsigned int prec) override
		// {
		// 	auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
		// 	val_pair.first.precision(prec);
		// }
				
				
		
				// /**
				// The computation of degree for Jacobians is challenging and not correctly implemented, so it is private.
				// */
				// int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
				// {
				//     return entry_node_->Degree(v);
				// }

	//       /**
	//       The computation of degree for Jacobians is challenging and not correctly implemented, so it is private.
	//       */
	//       int Degree(VariableGroup const& vars) const override
		// {
		// 	return entry_node_->Degree(vars);
		// }

		// /**
		//  Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	 //    */
		// std::vector<int> MultiDegree(VariableGroup const& vars) const override
		// {
			
		// 	std::vector<int> deg(vars.size());
		// 	for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
		// 	{
		// 		*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
		// 	}
		// 	return deg;
		// }

		// bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		// {
		// 	return entry_node_->IsHomogeneous(v);
		// }


		// /**
		// Check for homogeneity, with respect to a variable group.
		// */
		// bool IsHomogeneous(VariableGroup const& vars) const override
		// {
		// 	return entry_node_->IsHomogeneous(vars);
		// }

				std::tuple< std::pair<dbl,std::shared_ptr<Variable>>, std::pair<mpfr,std::shared_ptr<Variable>> > current_diff_variable_;
				Jacobian() = default;
		private:
				friend class boost::serialization::access;
		
				template <typename Archive>
				void serialize(Archive& ar, const unsigned version) {
						ar & boost::serialization::base_object<Function>(*this);
				}
		};
		
} // re: namespace node
} // re: namespace bertini



#endif
