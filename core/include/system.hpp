//This file is part of Bertini 2.0.
//
//system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
// system.hpp:  provides the bertini::system class.


#ifndef BERTINI_SYSTEM_H
#define BERTINI_SYSTEM_H

#include "mpfr_complex.hpp"

#include <vector>
#include "function_tree.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/number.hpp>

#include <assert.h>

#include <eigen3/Eigen/Dense>


namespace bertini {
	
	template<typename T>
	using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	
	/**
	 The fundamental polynomial system class for Bertini2.
	 */
	class System{
		
		
		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<Function>;
		using Var = std::shared_ptr<Variable>;
		
		
	public:
		
		void precision(unsigned new_precision)
		{
			for (auto iter : functions_) {
//				iter->precision(new_precision);
			}
			
			for (auto iter : subfunctions_) {
				//				iter->precision(new_precision);
			}
			
			for (auto iter : explicit_parameters_) {
				//				iter->precision(new_precision);
			}
			
			for (auto iter : variables_) {
				//				iter->precision(new_precision);
			}
			
			for (auto iter :implicit_parameters_) {
				//				iter->precision(new_precision);
			}
			
			for (auto iter : constants_) {
				//				iter->precision(new_precision);
				//				iter->eval<>()
			}
			
//			path_variable_->precision(new_precision);
			
			precision_ = new_precision;
		}
		
		
		
		
		/**
		 Get the number of functions in this system
		 */
		auto NumFunctions() const
		{
			return functions_.size();
		}
		
		
		
		/**
		 Evaluate the system.
		 */
		template<typename T>
		Vec<T> Eval()
		{
			Vec<T> value(NumFunctions()); // create vector with correct number of entries.
			
			{
				auto function_counter = 0;
				for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, function_counter++) {
					value(function_counter) = (*iter)->Eval<T>();
				}
			}
			
			return value;
		}
		
		
		/**
		 Set the values of the variables to be equal to the input values
		 */
		template<typename T>
		void SetVariables(Vec<T> new_values)
		{
			assert(new_values.size()== variables_.size());
			
			{
				auto counter = 0;
				for (auto iter=variables_.begin(); iter!=variables_.end(); iter++, counter++) {
					(*iter)->set_current_value(new_values(counter));
				}
			}
			
		}
		
		
		
		/**
		 Set the current value of the path variable.
		 */
		template<typename T>
		void SetPathVariable(T new_value)
		{
			path_variable_->set_current_value(new_value);
		}
		
		
		
		
		/**
		 For a system with implicitly defined parameters, set their values.  The values are determined externally to the system, and are tracked along with the variables.
		 */
		template<typename T>
		void SetImplicitParameters(Vec<T> new_values)
		{
			assert(new_values.size()== implicit_parameters_.size());
			
			{
				auto counter = 0;
				for (auto iter=implicit_parameters_.begin(); iter!=implicit_parameters_.end(); iter++, counter++) {
					(*iter)->set_current_value(new_values(counter));
				}
			}
		}
		
		
	private:
		
		
		std::vector< Fn > functions_;
		std::vector< Fn > subfunctions_;
		std::vector< Fn > explicit_parameters_;
		
		
		std::vector< Var > variables_;
		std::vector< Var > implicit_parameters_;
		
		
		Var path_variable_;
		
		std::vector< Fn > constants_;
		
		unsigned precision_;
		
		
		// i disagree with the inclusion of this, but the real necessity of it remains to be seen.  --dab
		Vec<bertini::complex> solutions_;
	};
	
}









#endif // for the ifndef include guards



