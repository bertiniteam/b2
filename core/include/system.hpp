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
		using Nd = std::shared_ptr<Node>;
		
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
			
			for (auto iter : constant_subfunctions_) {
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
		 Get the number of variables in this system
		 */
		auto NumVariables() const
		{
			return variables_.size();
		}
		
		/**
		 Get the number of constants in this system
		 */
		auto NumConstants() const
		{
			return constant_subfunctions_.size();
		}
		
		/**
		 Get the number of explicit parameters in this system
		 */
		auto NumParameters() const
		{
			return explicit_parameters_.size();
		}
		
		
		/**
		 Get the number of implicit parameters in this system
		 */
		auto NumImplicitParameters() const
		{
			return implicit_parameters_.size();
		}
		
		
		/**
		 Evaluate the system.
		 */
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values, const T & path_variable_value)
		{
			
			assert(variable_values.size()==NumVariables());
			
			// this function call traverses the entire tree, resetting everything.
			//
			// it has the unfortunate side effect of resetting constant functions, too.
			//
			// we need to work to correct this.
			for (auto iter : functions_) {
				iter->Reset();
			}
			
			SetVariables(variable_values);
			SetPathVariable(path_variable_value);
			
			
			Vec<T> value(NumFunctions()); // create vector with correct number of entries.
			
			{ // for scoping of the counter.
				auto counter = 0;
				for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, counter++) {
					value(counter) = (*iter)->Eval<T>();
				}
			}
			
			return value;
		}
		
		
		
		
		
		/**
		 Set the values of the variables to be equal to the input values
		 */
		template<typename T>
		void SetVariables(const Vec<T> & new_values)
		{
			assert(new_values.size()== variables_.size());
			
			{ // for scoping of the counter.
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
		
		
		
		
		
		
		
		
		/**
		 Add a variable group to the system.  The system will be homogenized with respect to this variable group, though this is not done at the time of this call.
		 */
		void AddVariableGroup(std::vector<Var> const& v)
		{
			variable_groups_.push_back(v);
			variables_.insert( variables_.end(), v.begin(), v.end() );
		}
		
		
		/**
		 Add a homogeneous (projective) variable group to the system.  The system must be homogeneous with respect to this group, though this is not verified at the time of this call.
		 */
		void AddHomVariableGroup(std::vector<Var> const& v)
		{
			hom_variable_groups_.push_back(v);
			variables_.insert( variables_.end(), v.begin(), v.end() );
		}
		
		
		/**
		 Add variables to the system which are in neither a regular variable group, nor in a homogeneous group.  This is likely used for user-defined systems, Classic userhomotopy: 1;.
		 */
		void AddUngroupedVariables(std::vector<Var> const& v)
		{
			ungrouped_variables_.insert( variables_.end(), v.begin(), v.end() );
			variables_.insert( variables_.end(), v.begin(), v.end() );
		}

		
		/**
		 Add an implicit parameter to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.
		 */
		void AddImplicitParameter(Var const& v)
		{
			implicit_parameters_.push_back(v);
		}
		
		
		/**
		 Add some implicit parameters to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.
		 */
		void AddImplicitParameters(std::vector<Var> const& v)
		{
			implicit_parameters_.insert( implicit_parameters_.end(), v.begin(), v.end() );
		}
		
		
		
		
		
		
		
		/**
		 Add an explicit parameter to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.
		 */
		void AddParameter(Fn const& F)
		{
			explicit_parameters_.push_back(F);
		}
		
		/**
		 Add some explicit parameters to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.
		 */
		void AddParameters(std::vector<Fn> const& v)
		{
			explicit_parameters_.insert( explicit_parameters_.end(), v.begin(), v.end() );
		}
		
		
		
		
		/**
		 Add a subfunction to the system.
		 */
		void AddSubfunction(Fn const& F)
		{
			subfunctions_.push_back(F);
		}
		
		/**
		 Add some subfunctions to the system.
		 */
		void AddSubfunctions(std::vector<Fn> const& v)
		{
			subfunctions_.insert( subfunctions_.end(), v.begin(), v.end() );
		}
		
		
		
		/**
		 Add a function to the system.
		 */
		void AddFunction(Fn const& F)
		{
			functions_.push_back(F);
		}
		
		/**
		 Add some functions to the system.
		 */
		void AddFunctions(std::vector<Fn> const& v)
		{
			functions_.insert( functions_.end(), v.begin(), v.end() );
		}
		
		
		
		
		
		/**
		 Add a constant function to the system.  Constants must not depend on anything which can vary -- they're constant!
		 */
		void AddConstant(Fn const& F)
		{
			constant_subfunctions_.push_back(F);
		}
		
		/**
		 Add some constant functions to the system.  Constants must not depend on anything which can vary -- they're constant!
		 */
		void AddConstants(std::vector<Fn> const& v)
		{
			constant_subfunctions_.insert( constant_subfunctions_.end(), v.begin(), v.end() );
		}
		
		
		
		
		
		/**
		 Add a variable as the Path Variable to a System.  Will overwrite any previously declared path variable.
		 */
		void AddPathVariable(Var const& v)
		{
			path_variable_ = v;
		}
		
		
		
		
		/**
		 Overloaded operator for printing to an arbirtary out stream.
		 */
		friend std::ostream& operator <<(std::ostream& out, const System & s)
		{
			out << "system:\n\n";
			out << s.NumVariables() << " variables:\n";
			for (auto iter : s.variables_) {
				out << (iter)->name() << "\n";
			}
			out << "\n";
			
			out << s.NumFunctions() << " functions:\n";
			for (auto iter : s.functions_) {
				out << (iter)->name() << "\n";
				out << *iter << "\n";
			}
			out << "\n";
			
			
			if (s.NumParameters()) {
				out << s.NumParameters() << " explicit parameters:\n";
				for (auto iter : s.explicit_parameters_) {
					out << (iter)->name() << "\n";
					out << *iter << "\n";
				}
				out << "\n";
			}
			
			
			if (s.NumConstants()) {
				out << s.NumConstants() << " constants:\n";
				for (auto iter : s.constant_subfunctions_) {
					out << (iter)->name() << "\n";
					out << *iter << "\n";
				}
				out << "\n";
			}
			
			if (s.path_variable_) {
				out << "path variable defined.  named " << s.path_variable_->name() << "\n";
			}
			
			return out;
		}
		
		
	private:
		
		
		std::vector<Var> ungrouped_variables_;
		std::vector<std::vector<Var> > variable_groups_;
		std::vector<std::vector<Var> > hom_variable_groups_;
		
		
		std::vector< Fn > functions_;
		std::vector< Fn > subfunctions_;
		std::vector< Fn > explicit_parameters_;
		
		
		std::vector< Var > variables_;
		std::vector< Var > implicit_parameters_;
		
		
		Var path_variable_;
		
		std::vector< Fn > constant_subfunctions_;
		
		unsigned precision_;
		
		
		// i disagree with the inclusion of this, but the real necessity of it remains to be seen.  --dab
		Vec<bertini::complex> solutions_;
	};
	
}









#endif // for the ifndef include guards



