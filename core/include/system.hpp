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

	
	/**
	 The fundamental polynomial system class for Bertini2.
	 */
	class System{
		

	public:

		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<Function>;
		using Var = std::shared_ptr<Variable>;
		using Nd = std::shared_ptr<Node>;
		using Jac = std::shared_ptr<Jacobian>;
		
		System() : is_differentiated_(false), have_path_variable_(false)
		{}


		void precision(unsigned new_precision);




		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values);

		/**
		 Evaluate the system.
		 */
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values, const T & path_variable_value);



		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values);



		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values, const T & path_variable_value);

	

		void Homogenize();


		bool IsHomogeneous() const;


		bool IsPolynomial() const;


		//////////////////
		//
		//  Nummers --   functions which get the numbers of things.
		//
		//////////////////


		/**
		 Get the number of functions in this system
		 */
		size_t NumFunctions() const;

		/**
		 Get the number of variables in this system
		 */
		size_t NumVariables() const;

		/**
		 Get the number of variable groups in the system
		*/
		 size_t NumVariableGroups() const;

		/**
		 Get the number of homogeneous variable groups in the system
		*/
		 size_t NumHomVariableGroups() const;


		/**
		 Get the number of constants in this system
		 */
		size_t NumConstants() const;

		/**
		 Get the number of explicit parameters in this system
		 */
		size_t NumParameters() const;


		/**
		 Get the number of implicit parameters in this system
		 */
		size_t NumImplicitParameters() const;


		


		///////////////////
		//
		// Setters -- templated.
		//
		///////////////

		/**
		 Set the values of the variables to be equal to the input values
		 */
		template<typename T>
		void SetVariables(const Vec<T> & new_values);



		/**
		 Set the current value of the path variable.
		 */
		template<typename T>
		void SetPathVariable(T new_value);


		


		/**
		 For a system with implicitly defined parameters, set their values.  The values are determined externally to the system, and are tracked along with the variables.
		 */
		template<typename T>
		void SetImplicitParameters(Vec<T> new_values);








		//////////////////
		//
		//  Adders  --   functions which add things to the system.
		//
		//////////////////


		/**
		 Add a variable group to the system.  The system will be homogenized with respect to this variable group, though this is not done at the time of this call.
		 */
		void AddVariableGroup(VariableGroup const& v);


		/**
		 Add a homogeneous (projective) variable group to the system.  The system must be homogeneous with respect to this group, though this is not verified at the time of this call.
		 */
		void AddHomVariableGroup(VariableGroup const& v);



		/**
		 Add variables to the system which are in neither a regular variable group, nor in a homogeneous group.  This is likely used for user-defined systems, Classic userhomotopy: 1;.
		 */
		void AddUngroupedVariable(Var const& v);


		/**
		 Add variables to the system which are in neither a regular variable group, nor in a homogeneous group.  This is likely used for user-defined systems, Classic userhomotopy: 1;.
		 */
		void AddUngroupedVariables(VariableGroup const& v);


		/**
		 Add an implicit parameter to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.
		 */
		void AddImplicitParameter(Var const& v);


		/**
		 Add some implicit parameters to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.
		 */
		void AddImplicitParameters(VariableGroup const& v);




		/**
		 Add an explicit parameter to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.
		 */
		void AddParameter(Fn const& F);

		/**
		 Add some explicit parameters to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.
		 */
		void AddParameters(std::vector<Fn> const& v);



		/**
		 Add a subfunction to the system.
		 */
		void AddSubfunction(Fn const& F);

		/**
		 Add some subfunctions to the system.
		 */
		void AddSubfunctions(std::vector<Fn> const& v);



		/**
		 Add a function to the system.
		 */
		void AddFunction(Fn const& F);

		/**
		 Add a function to the system.
		 */
		void AddFunction(Nd const& N);


		/**
		 Add some functions to the system.
		 */
		void AddFunctions(std::vector<Fn> const& v);







		/**
		 Add a constant function to the system.  Constants must not depend on anything which can vary -- they're constant!
		 */
		void AddConstant(Fn const& F);


		/**
		 Add some constant functions to the system.  Constants must not depend on anything which can vary -- they're constant!
		 */
		void AddConstants(std::vector<Fn> const& v);




		/**
		 Add a variable as the Path Variable to a System.  Will overwrite any previously declared path variable.
		 */
		void AddPathVariable(Var const& v);


		bool HavePathVariable() const;








		/////////////// TESTING ////////////////////
		/**
		 Get a function by its index.  This is just as scary as you think it is.  It is up to you to make sure the function at this index exists.
		*/
        auto function(unsigned index = 0) const
        {
            return functions_[index];
        }

        /**
		 Get the variables in the problem.
        */
        auto variables() const
        {
            return variables_;
        }

        /**
		 Get the variable groups in the problem.
        */
        auto variableGroups() const
        {
            return variable_groups_;
        }

		/////////////// TESTING ////////////////////




		/**
		 Get the degrees of the functions in the system, with respect to all variables.
		*/
		 std::vector<int> Degrees() const;

		 /**
		 Get the degrees of the functions in the system, with respect to a group of variables.
		*/
		 std::vector<int> Degrees(VariableGroup const& vars) const;

		/**
		 Sort the functions so they are in decreasing order by degree
		*/
		void ReorderFunctionsByDegreeDecreasing();


		/**
		 Sort the functions so they are in decreasing order by degree
		*/
		void ReorderFunctionsByDegreeIncreasing();

		/**
		 Overloaded operator for printing to an arbirtary out stream.
		 */
		friend std::ostream& operator <<(std::ostream& out, const System & s);


		/**
		 Clear the entire structure of variables in a system.
		*/
		void ClearVariables();


		/**
		 Copy the entire structure of variables from within one system to another.
		  This copies everything -- ungrouped variables, variable groups, homogenizing variables, the path variable, the ordering of the variables.s
		*/ 
		void CopyVariableStructure(System const& other);
        

		System operator+=(System const& rhs);
		friend System operator+(System lhs, System const& rhs);

		System operator*=(std::shared_ptr<Node> const& N);
		friend System operator*(System s, std::shared_ptr<Node> const&  N);
		friend System operator*(std::shared_ptr<Node> const&  N, System const& s);
	private:

		std::vector<Var> ungrouped_variables_;
		std::vector< VariableGroup > variable_groups_;
		std::vector< VariableGroup > hom_variable_groups_;

		VariableGroup homogenizing_variables_;

		std::vector< Fn > functions_;
		std::vector< Fn > subfunctions_;
		std::vector< Fn > explicit_parameters_;


		VariableGroup variables_;
		VariableGroup implicit_parameters_;


		Var path_variable_;
		bool have_path_variable_;


		std::vector< Fn > constant_subfunctions_;

		std::vector< Jac > jacobian_;
		bool is_differentiated_;


		unsigned precision_;

	};




	




}









#endif // for the ifndef include guards



