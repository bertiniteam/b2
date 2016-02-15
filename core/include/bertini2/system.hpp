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

#include "bertini2/mpfr_complex.hpp"

#include <vector>
#include "bertini2/function_tree.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/number.hpp>

#include <assert.h>

#include <Eigen/Dense>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>


namespace bertini {

	
	/**
	\brief The fundamental polynomial system class for Bertini2.
	
	 The fundamental polynomial system class for Bertini2.

	 Other System types are derived from this, but this class is not abstract.
	 */
	class System{
		

	public:

		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<node::Function>;
		using Var = std::shared_ptr<node::Variable>;
		using Nd = std::shared_ptr<node::Node>;
		using Jac = std::shared_ptr<node::Jacobian>;
		
		/**
		The default constructor for a system
		*/
		System() : is_differentiated_(false), have_path_variable_(false)
		{}


		/**
		Change the precision of the entire system's functions, subfunctions, and all other nodes.

		\param new_precision The new precision, in digits, to work in.  This only affects the mpfr types, not double.  To use low-precision (doubles), use that number type in the templated functions.
		*/
		void precision(unsigned new_precision);


		/**
		 \brief Compute and internally store the symbolic Jacobian of the system.
		*/
		void Differentiate();

		/**
		Evaluate the system, provided the system has no path variable defined.

		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values)
		{

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (have_path_variable_)
				throw std::runtime_error("not using a time value for evaluation of system, but path variable IS defined.");


			// this function call traverses the entire tree, resetting everything.
			//
			// TODO: it has the unfortunate side effect of resetting constant functions, too.
			//
			// we need to work to correct this.
			for (auto iter : functions_) {
				iter->Reset();
			}

			SetVariables(variable_values);


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
		 Evaluate the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.
		 */
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values, const T & path_variable_value)
		{

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");


			// this function call traverses the entire tree, resetting everything.
			//
			// TODO: it has the unfortunate side effect of resetting constant functions, too.
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
		Evaluate the Jacobian matrix of the system, provided the system has no path variable defined.

		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values)
		{
			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (have_path_variable_)
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");


			if (!is_differentiated_)
				Differentiate();


			SetVariables(variable_values);

			auto vars = Variables(); //TODO: replace this with something that peeks directly into the variables without this copy.

			Mat<T> J(NumFunctions(), NumVariables());
			for (int ii = 0; ii < NumFunctions(); ++ii)
				for (int jj = 0; jj < NumVariables(); ++jj)
					J(ii,jj) = jacobian_[ii]->EvalJ<T>(vars[jj]);

			return J;
		}


		/**
		 Evaluate the Jacobian of the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \return The Jacobian matrix.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 */
		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values, const T & path_variable_value)
		{
			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");


			if (!is_differentiated_)
				Differentiate();

			SetVariables(variable_values);
			SetPathVariable(path_variable_value);

			auto vars = Variables(); //TODO: replace this with something that peeks directly into the variables without this copy.
			
			Mat<T> J(NumFunctions(), NumVariables());
			for (int ii = 0; ii < NumFunctions(); ++ii)
				for (int jj = 0; jj < NumVariables(); ++jj)
					J(ii,jj) = jacobian_[ii]->EvalJ<T>(vars[jj]);

			return J;

		}

	
		/**
		Homogenize the system, adding new homogenizing variables for each VariableGroup defined for the system.

		\throws std::runtime_error, if the system is not polynomial, has a mismatch on the number of homogenizing variables and the number of variable groups (this would result from a partially homogenized system), or the homogenizing variable names somehow get screwed up by having duplicates.
		*/
		void Homogenize();

		/**
		Checks whether a system is homogeneous, overall.  This means with respect to each variable group (including homogenizing variable if defined), homogeneous variable group, and ungrouped variables, if defined.

		\throws std::runtime_error, if the number of homogenizing variables does not match the number of variable_groups.
		\return true if homogeneous, false if not
		*/
		bool IsHomogeneous() const;

		/**
		Checks whether a system is polynomial, overall.  This means with respect to each variable group (including homogenizing variable if defined), homogeneous variable group, and ungrouped variables, if defined.

		\throws std::runtime_error, if the number of homogenizing variables does not match the number of variable_groups.
		\return true if polynomial, false if not.
		*/
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
		 Get the number of *homogenizing* variables in this system
		 */
		size_t NumHomVariables() const;
		/**
		 Get the number of variable groups in the system
		*/
		 size_t NumVariableGroups() const;


		/**
		 get the number of variables which are ungrouped.
		 */
		 size_t NumUngroupedVariables() const;
		 
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

		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 \throws std::runtime_error if the number of variables doesn't match.
		 The ordering of the variables matters.  The standard ordering is 1) variable groups, with homogenizing variable first. 2) homogeneous variable groups. 3) ungrouped variables.

		 The path variable is not considered a variable for this operation.  It is set separately.
		 
		 \param new_values The new updated values for the variables.

		 \see SetPathVariable
		 \see Variables
		 */
		template<typename T>
		void SetVariables(const Vec<T> & new_values)
		{
			assert(new_values.size()== NumVariables());

			auto vars = Variables();

			auto counter = 0;

			for (auto iter=vars.begin(); iter!=vars.end(); iter++, counter++) {
				(*iter)->set_current_value(new_values(counter));
			}
		}



		/**
		 Set the current value of the path variable.
		
		 \throws std::runtime_error, if a path variable is not defined.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

		 \param new_value The new updated values for the path variable.
		 */
		template<typename T>
		void SetPathVariable(T new_value)
		{
			if (!have_path_variable_)
				throw std::runtime_error("trying to set the value of the path variable, but one is not defined for this system");

			path_variable_->set_current_value(new_value);
		}


		


		/**
		 For a system with implicitly defined parameters, set their values.  The values are determined externally to the system, and are tracked along with the variables.
		 
		 \param new_values The new updated values for the implicit parameters.
		 */
		template<typename T>
		void SetImplicitParameters(Vec<T> new_values)
		{
			if (new_values.size()!= implicit_parameters_.size())
				throw std::runtime_error("trying to set implicit parameter values, but there is a size mismatch");

			size_t counter = 0;
			for (auto iter=implicit_parameters_.begin(); iter!=implicit_parameters_.end(); iter++, counter++)
				(*iter)->set_current_value(new_values(counter));

		}








		//////////////////
		//
		//  Adders  --   functions which add things to the system.
		//
		//////////////////


		/**
		 Add a variable group to the system.  The system may be homogenized with respect to this variable group, though this is not done at the time of this call.
		 
		 \param v The variable group to add.
		 */
		void AddVariableGroup(VariableGroup const& v);


		/**
		 Add a homogeneous (projective) variable group to the system.  The system must be homogeneous with respect to this group, though this is not verified at the time of this call.
		 
		 \param v The variable group to add.
		 */
		void AddHomVariableGroup(VariableGroup const& v);



		/**
		 Add variables to the system which are in neither a regular variable group, nor in a homogeneous group.  This is likely used for user-defined systems, Classic userhomotopy: 1;.
		 
		 \param v The variable to add.
		 */
		void AddUngroupedVariable(Var const& v);


		/**
		 Add variables to the system which are in neither a regular variable group, nor in a homogeneous group.  This is likely used for user-defined systems, Classic userhomotopy: 1;.
		 
		 \param v The variables to add.
		 */
		void AddUngroupedVariables(VariableGroup const& v);


		/**
		 Add an implicit parameter to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.

		 \param v The implicit parameter to add.
		 */
		void AddImplicitParameter(Var const& v);


		/**
		 Add some implicit parameters to the system.  Implicit parameters are advanced by the tracker akin to variable advancement.

		 \param v The implicit parameters to add.
		 */
		void AddImplicitParameters(VariableGroup const& v);




		/**
		 Add an explicit parameter to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.

		 \param F The parameter to add.
		 */
		void AddParameter(Fn const& F);

		/**
		 Add some explicit parameters to the system.  Explicit parameters should depend only on the path variable, though this is not checked in this function.

		 \param F The parameters to add.
		 */
		void AddParameters(std::vector<Fn> const& F);



		/**
		 Add a subfunction to the system.

		 Tacks them onto the end of the system.

		 \param F The subfunction to add.
		 */
		void AddSubfunction(Fn const& F);

		/**
		 Add some subfunctions to the system.

		 Tacks them onto the end of the system.

		 \param F The subfunctions to add.
		 */
		void AddSubfunctions(std::vector<Fn> const& F);



		/**
		 Add a function to the system.

		 \param F The function to add.
		 */
		void AddFunction(Fn const& F);

		/**
		 Add a function to the system.

		 \param F The function to add.
		 */
		void AddFunction(Nd const& F);


		/**
		 Add some functions to the system.

		 \param F The functions to add.
		 */
		void AddFunctions(std::vector<Fn> const& F);







		/**
		 Add a constant function to the system.  Constants must not depend on anything which can vary -- they're constant!

		 \param C The constant to add.
		 */
		void AddConstant(Fn const& C);


		/**
		 Add some constant functions to the system.  Constants must not depend on anything which can vary -- they're constant!

		 \param C The constants to add.
		 */
		void AddConstants(std::vector<Fn> const& C);




		/**
		 Add a variable as the Path Variable to a System.  Will overwrite any previously declared path variable.

		 \param v The new path variable.
		 */
		void AddPathVariable(Var const& v);


		/**
		Query whether a path variable is set for this system

		\return true if have a path variable, false if not.
		*/
		bool HavePathVariable() const;






		/**
		 Get the variables in the problem.

		 This function returns the variables in standard ordering.

		 1) variable groups first, lead by their respective homogenizing variable, if defined.
		 2) hom_variable_groups second.
		 3) ungrouped variables

		 The order in which variables and their groups are added to a system impacts this ordering.  Bertini will ensure internal consistency.  It is up to the user to make sure the ordering is what they want for any post-processing or interpretations of results.

		 \throws std::runtime_error, if there is a mismatch between the number of homogenizing variables and the number of variable_groups.  This would happen if a system is homogenized, and then more stuff is added to it.  
		*/
		VariableGroup Variables() const;




		/////////////// TESTING ////////////////////
		/**
		 Get a function by its index.  This is just as scary as you think it is.  It is up to you to make sure the function at this index exists.
		*/
		auto Function(unsigned index) const
		{
			return functions_[index];
		}

		

		/**
		 Get the variable groups in the problem.
		*/
		auto VariableGroups() const
		{
			return variable_groups_;
		}


		/**
		 Get the homogeneous variable groups in the problem.
		*/
		auto HomVariableGroups() const
		{
			return hom_variable_groups_;
		}
		/////////////// TESTING ////////////////////




		/**
		 Get the degrees of the functions in the system, with respect to all variables.

		 \return A vector containing the degrees of the functions.  Negative numbers indicate the function is non-polynomial.
		*/
		 std::vector<int> Degrees() const;

		 /**
		 Get the degrees of the functions in the system, with respect to a group of variables.

		 \return A vector containing the degrees of the functions.  Negative numbers indicate the function is non-polynomial.
		 \param vars A group of variables with respect to which you wish to compute degrees.  Needs not be a group with respect to the system.
		*/
		 std::vector<int> Degrees(VariableGroup const& vars) const;

		/**
		 Sort the functions so they are in DEcreasing order by degree
		*/
		void ReorderFunctionsByDegreeDecreasing();


		/**
		 Sort the functions so they are in INcreasing order by degree
		*/
		void ReorderFunctionsByDegreeIncreasing();

		/**
		 Overloaded operator for printing to an arbirtary out stream.
		 */
		friend std::ostream& operator <<(std::ostream& out, const System & s);


		/**
		 Clear the entire structure of variables in a system.  Reconstructing it is up to you.
		*/
		void ClearVariables();


		/**
		 Copy the entire structure of variables from within one system to another.
		  This copies everything -- ungrouped variables, variable groups, homogenizing variables, the path variable, the ordering of the variables.

		  \param other Another system from which to copy the variable structure.  

		  This operation does NOT affect the functions in any way.  It is up to the user to make the functions actually depend on these variables.
		*/ 
		void CopyVariableStructure(System const& other);
		

		/**
		Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.

		*/
		System operator+=(System const& rhs);

		/**
		Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.
		
		*/
		friend System operator+(System lhs, System const& rhs);

		/**
		Multiply a system by an arbitrary node.  Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		System operator*=(Nd const& N);

		/**
		Multiply a system by an arbitrary node.  Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend System operator*(System s, Nd const&  N);

		/**
		Multiply a system by an arbitrary node.  Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend System operator*(Nd const&  N, System const& s);
	private:

		VariableGroup ungrouped_variables_; ///< ungrouped variable nodes.  Not in an affine variable group, not in a projective group.  Just hanging out, being a variable.
		std::vector< VariableGroup > variable_groups_; ///< Affine variable groups.  When system is homogenized, will have a corresponding homogenizing variable.
		std::vector< VariableGroup > hom_variable_groups_; ///< Homogeneous or projective variable groups.  System SHOULD be homogeneous with respect to these.  

		VariableGroup homogenizing_variables_; ///< homogenizing variables for the variable_groups.  


		bool have_path_variable_; ///< Whether we have the variable or not.
		Var path_variable_; ///< the single path variable for this system.  Sometimes called time.
		
		VariableGroup implicit_parameters_; ///< Implicit parameters.  These don't depend on anything, and will be moved from one parameter point to another by the tracker.  They should be algebraically constrained by some equations.
		std::vector< Fn > explicit_parameters_; ///< Explicit parameters.  These should be functions of the path variable only, NOT of other variables.  

		std::vector< Fn > constant_subfunctions_; ///< degree-0 functions, depending on neither variables nor the path variable.
		std::vector< Fn > subfunctions_; ///< Any declared subfunctions for the system.  Can use these to ensure that complicated repeated structures are only created and evaluated once.
		std::vector< Fn > functions_; ///< The system's functions.
		
		
		std::vector< Jac > jacobian_; ///< The generated functions from differentiation.  Created when first call for a Jacobian matrix evaluation.
		bool is_differentiated_; ///< indicator for whether the jacobian tree has been populated.


		unsigned precision_; ///< the current working precision of the system 


		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & ungrouped_variables_;
			ar & variable_groups_;
			ar & hom_variable_groups_;
			ar & homogenizing_variables_;
			
			ar & have_path_variable_;
			ar & path_variable_;

			ar & implicit_parameters_;
			ar & explicit_parameters_;

			ar & constant_subfunctions_;
			ar & subfunctions_;
			ar & functions_;

			ar & is_differentiated_;
			ar & jacobian_;
		}

	};




	




}









#endif // for the ifndef include guards



