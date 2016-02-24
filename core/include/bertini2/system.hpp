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


/**
\file system.hpp 

\brief Provides the bertini::system class.
*/

#ifndef BERTINI_SYSTEM_HPP
#define BERTINI_SYSTEM_HPP

#include "mpfr_complex.hpp"
#include "mpfr_extensions.hpp"
#include "eigen_extensions.hpp"

#include <vector>
#include "function_tree.hpp"

#include "patch.hpp"

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/number.hpp>

#include <assert.h>



#include <boost/type_index.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>

#include "limbo.hpp"
#include <Eigen/Dense>

namespace bertini {

	

	/**
	\brief The fundamental polynomial system class for Bertini2.
	
	 The fundamental polynomial system class for Bertini2.

	 Other System types are derived from this, but this class is not abstract.
	 */
	class System{
		
	private:
		enum class OrderingChoice{FIFO, AffHomUng};

	public:
		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<node::Function>;
		using Var = std::shared_ptr<node::Variable>;
		using Nd = std::shared_ptr<node::Node>;
		using Jac = std::shared_ptr<node::Jacobian>;
		
		/**
		\brief The default constructor for a system.
		*/
		System() : is_differentiated_(false), have_path_variable_(false), have_ordering_(false), ordering_(OrderingChoice::FIFO), precision_(mpfr_float::default_precision()), is_patched_(false)
		{}

		/** 
		\brief The copy operator, creates a system from a string using the Bertini parser for Bertini classic syntax.
		*/
		System(std::string const& input);
		
		/** 
		\brief The copy operator
		*/
		System(System const& other);

		/** 
		\brief The move copy operator
		*/
		System(System && other) : System()
		{
			swap(*this, other);
		}

		/** 
		\brief The assignment operator
		*/
		System& operator=(System other);

		/**
		\brief The move assignment operator
		*/
		System& operator=(System && other) = default;

		/**
		The free swap function for systems.
		*/
		friend void swap(System & a, System & b);

		/**
		Change the precision of the entire system's functions, subfunctions, and all other nodes.

		\param new_precision The new precision, in digits, to work in.  This only affects the mpfr types, not double.  To use low-precision (doubles), use that number type in the templated functions.
		*/
		void precision(unsigned new_precision) const;

		/**
		\brief Get the current precision of a system.
		*/
		unsigned precision() const
		{
			return precision_;
		}

		/**
		 \brief Compute and internally store the symbolic Jacobian of the system.
		*/
		void Differentiate() const;


		/**
		\brief Evaluate the system using the previously set variable (and time) values.  

		It is up to YOU to ensure that the system's variables (and path variable) has been set prior to this function call.

		\return The function values of the system
		*/ 
		template<typename T>
		Vec<T> Eval() const
		{


			// the Reset() function call traverses the entire tree, resetting everything.
			// TODO: it has the unfortunate side effect of resetting constant functions, too.
			for (auto iter : functions_) 
				iter->Reset();

			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.

			unsigned counter(0);
			for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, counter++) {
				function_values(counter) = (*iter)->Eval<T>();
			}

			if (IsPatched())
				patch_.Eval(function_values,std::get<Vec<T> >(current_variable_values_));// .segment(NumFunctions(),NumTotalVariableGroups())

			return function_values;
		}


		/**
		\brief Evaluate the system, provided the system has no path variable defined.

		Causes the current variable values to be set in the system.  Resets the function tree's stored numbers.  


		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values) const
		{

			if (variable_values.size()!=NumVariables())
			{
				std::stringstream ss;
				ss << "trying to evaluate system, but number of input variables (" << variable_values.size() << ") doesn't match number of system variables (" << NumVariables() << ").";
				throw std::runtime_error(ss.str());
			}
			if (have_path_variable_)
				throw std::runtime_error("not using a time value for evaluation of system, but path variable IS defined.");

			SetVariables(variable_values);

			return Eval<T>();
		}



		/**
		 Evaluate the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \todo The Eval() function for systems has the unfortunate side effect of resetting constant functions.  Modify the System class so that only certain parts of the tree get reset.
		 */
		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values, const T & path_variable_value) const
		{

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");

			SetVariables(variable_values);
			SetPathVariable(path_variable_value);

			return Eval<T>();
		}


		/**
		Evaluate the Jacobian matrix of the system, using the previous space and time values.

		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		*/
		template<typename T>
		Mat<T> Jacobian() const
		{
			auto vars = Variables(); //TODO: replace this with something that peeks directly into the variables without this copy.

			if (!is_differentiated_)
				Differentiate();

			Mat<T> J(NumTotalFunctions(), NumVariables());
			for (int ii = 0; ii < NumFunctions(); ++ii)
				for (int jj = 0; jj < NumVariables(); ++jj)
					J(ii,jj) = jacobian_[ii]->EvalJ<T>(vars[jj]);
			if (IsPatched())
				patch_.Jacobian(J,std::get<Vec<T> >(current_variable_values_));

			return J;
		}

		/**
		Evaluate the Jacobian matrix of the system, provided the system has no path variable defined.

		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values) const
		{
			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (have_path_variable_)
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");

			SetVariables(variable_values);

			return Jacobian<T>();
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
		Mat<T> Jacobian(const Vec<T> & variable_values, const T & path_variable_value) const
		{
			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");

			SetVariables(variable_values);
			SetPathVariable(path_variable_value);

			return Jacobian<T>();
		}

		
		/**
		\brief Compute the time-derivative is a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename T>
		Vec<T> TimeDerivative(const Vec<T> & variable_values, const T & path_variable_value) const
		{
			if (!HavePathVariable())
				throw std::runtime_error("computing time derivative of system with no path variable defined");

			if (!is_differentiated_)
				Differentiate();

			SetVariables(variable_values);
			SetPathVariable(path_variable_value);

			Vec<T> ds_dt(NumTotalFunctions());
			
			for (int ii = 0; ii < NumFunctions(); ++ii)
				ds_dt(ii) = jacobian_[ii]->EvalJ<T>(path_variable_);		

			if (IsPatched())
				for (int ii = 0; ii < NumTotalVariableGroups(); ++ii)
					ds_dt(ii+NumFunctions()) = T(0);

			return ds_dt;
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
		 Get the number of functions in this system, excluding patches.
		 */
		size_t NumFunctions() const;

		/**
		Get the number of patches in this system.
		*/
		size_t NumPatches() const
		{
			return patch_.NumVariableGroups();
		}

		/**
		Get the total number of functions, including patches
		*/
		size_t NumTotalFunctions() const;

		/**
		 Get the total number of variables in this system, including homogenizing variables.
		 */
		size_t NumVariables() const;

		/**
		 Get the number of variables in this system, NOT including homogenizing variables.
		*/
		size_t NumNaturalVariables() const;

		/**
		 Get the number of *homogenizing* variables in this system
		 */
		size_t NumHomVariables() const;

		/**
		Get the total number of variable groups in the system, including both affine and homogenous.  Ignores the ungrouped variables, because they are not in any group.
		*/
		size_t NumTotalVariableGroups() const;

		/**
		 Get the number of affine variable groups in the system
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

		 The ordering of the variables matters.  

		 * The AffHomUng ordering is 1) variable groups, with homogenizing variable first. 2) homogeneous variable groups. 3) ungrouped variables.
		 * The FIFO ordering uses the order in which the variable groups were added.

		 The path variable is not considered a variable for this operation.  It is set separately.
		 
		 \param new_values The new updated values for the variables.

		 \see SetPathVariable
		 \see Variables
		 */
		template<typename T>
		void SetVariables(const Vec<T> & new_values) const
		{
			if (new_values.size()!= NumVariables())
				throw std::runtime_error("variable vector of different length from system-owned variables in SetVariables");

			#ifndef BERTINI_DISABLE_ASSERTS
			if (Precision(new_values(0))!=DoublePrecision())
				assert(Precision(new_values(0)) == precision() && "precision of input point in SetVariables must match the precision of the system.");
			#endif

			auto vars = Variables();

			auto counter = 0;

			for (auto iter=vars.begin(); iter!=vars.end(); iter++, counter++) {
				(*iter)->set_current_value(new_values(counter));
			}

			std::get<Vec<T> >(current_variable_values_) = new_values;
		}



		/**
		 Set the current value of the path variable.
		
		 \throws std::runtime_error, if a path variable is not defined.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

		 \param new_value The new updated values for the path variable.
		 */
		template<typename T>
		void SetPathVariable(T new_value) const
		{
			if (!have_path_variable_)
				throw std::runtime_error("trying to set the value of the path variable, but one is not defined for this system");

			path_variable_->set_current_value(new_value);
		}


		


		/**
		 For a system with implicitly defined parameters, set their values.  The values are determined externally to the system, and are tracked along with the variables.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 \param new_values The new updated values for the implicit parameters.
		 */
		template<typename T>
		void SetImplicitParameters(Vec<T> new_values) const
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
		/brief Query the current variable ordering.
		*/
		OrderingChoice Ordering() const
		{
			return ordering_;
		}

		/**
		 Order the variables in the problem, 1 affine, 2 homogeneous, 3 ungrouped.

		 This function returns the variables in standard ordering.

		 1) variable groups first, lead by their respective homogenizing variable, if defined.
		 2) hom_variable_groups second.
		 3) ungrouped variables

		 \throws std::runtime_error, if there is a mismatch between the number of homogenizing variables and the number of variable_groups.  This would happen if a system is homogenized, and then more stuff is added to it.  
        */
        VariableGroup AffHomUngVariableOrdering() const;



        /**
		 Order the variables, by the order in which the groups were added.

		 This function returns the variables in First In First Out (FIFO) ordering.

		 Homogenizing variables precede affine variable groups, so that groups always are grouped together.

		 \throws std::runtime_error, if there is a mismatch between the number of homogenizing variables and the number of variable_groups.  This would happen if a system is homogenized, and then more stuff is added to it.  
        */
        VariableGroup FIFOVariableOrdering() const;


        /**
		\brief Choose to use a canonical ordering for variables.  

		The system must be homogenized correctly for this to succeed.  Otherwise the underlying call to FIFOVariableOrdering may throw.

		\see AffHomUngariableOrdering
        */
        void SetAffHomUngVariableOrdering();


        /**
		\brief Choose to use order-of-addition ordering for variables.  

		The system must be homogenized correctly for this to succeed.  Otherwise the underlying call to FIFOVariableOrdering may throw.

		\see FIFOVariableOrdering
        */
        void SetFIFOVariableOrdering();

        /**
		 Get the variables in the problem; they must have already been ordered.

		 \see SetAffHomUngVariableOrdering
		 \see SetFIFOVariableOrdering
		*/
		VariableGroup Variables() const;

		/**
		\brief Get an affine variable group the class has defined.

		It is up to you to ensure this group exists.
		*/
		VariableGroup const& AffineVariableGroup(size_t index) const
		{
			return variable_groups_[index];
		}
		/**
		\brief Get the sizes of the variable groups, according to the current ordering
		*/
		std::vector<unsigned> VariableGroupSizes() const
		{
			if (Ordering()==OrderingChoice::FIFO)
				return VariableGroupSizesFIFO();
			else
				return VariableGroupSizesAffHomUng();
		}

        /**
		\brief Dehomogenize a point, using the variable grouping / structure of the system.
		
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

		\throws std::runtime_error, if there is a mismatch between the number of variables in the input point, and the total number of var
        */
        template<typename T>
	    Vec<T> DehomogenizePoint(Vec<T> const& x) const
	        {

	        	if (x.size()!=NumVariables())
	        		throw std::runtime_error("dehomogenizing point with incorrect number of coordinates");

	        	if (!have_ordering_)
	    			ConstructOrdering();



	    		switch (Ordering())
	    		{
	    			case OrderingChoice::AffHomUng:
	    			{
	    				return DehomogenizePointAffHomUngOrdering(x);
	    				break;
	    			}
	    			case OrderingChoice::FIFO:
	    			{
	    				return DehomogenizePointFIFOOrdering(x);
	    				break;
	    			}
	    			default: //ordering_==OrderingChoice::???
	    			{
	    				throw std::runtime_error("invalid ordering choice");
	    				break;
	    			}
	    		}
	        }



		/////////////// TESTING ////////////////////
		/**
		 Get a function by its index.  

		 This is just as scary as you think it is.  It is up to you to make sure the function at this index exists.
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
		Compute an estimate of an upper bound of the absolute values of the coefficients in the system.
		
		\param num_evaluations The number of times to compute this estimate.  Default is 1.
		\returns An upper bound on the absolute values of the coefficients.
		*/
        mpfr_float CoefficientBound(unsigned num_evaluations=1) const;


        /**
         \brief Compute an upper bound on the degree of the system.  

         This number will be wrong if the system is non-polynomial, because degree for non-polynomial systems is not defined.
         */
        int DegreeBound() const;

		/**
		 \brief Get the degrees of the functions in the system, with respect to all variables.

		 \return A vector containing the degrees of the functions.  Negative numbers indicate the function is non-polynomial.
		*/
		 std::vector<int> Degrees() const;

		 /**
		 \brief Get the degrees of the functions in the system, with respect to a group of variables.

		 \return A vector containing the degrees of the functions.  Negative numbers indicate the function is non-polynomial.
		 \param vars A group of variables with respect to which you wish to compute degrees.  Needs not be a group with respect to the system.
		*/
		 std::vector<int> Degrees(VariableGroup const& vars) const;

		/**
		 \brief Sort the functions so they are in DEcreasing order by degree
		*/
		void ReorderFunctionsByDegreeDecreasing();


		/**
		 \brief Sort the functions so they are in INcreasing order by degree
		*/
		void ReorderFunctionsByDegreeIncreasing();





		/////////////
		//
		//  Functions regarding patches.
		//
		//////////////





		/**
		\brief Let the system patch itself by the selected variable ordering, using the current variable groups.  

		Homogeneous variable groups and affine variable groups will be supplied a patch equation.  Ungrouped variables will not.

		\todo Add example code for how to use this function.
		*/
		void AutoPatch()
		{
			if (Ordering()==OrderingChoice::FIFO)
				AutoPatchFIFO();
			else if (Ordering()==OrderingChoice::AffHomUng)
				AutoPatchAffHomUng();
		}

		/**
		\brief Copy the patches from another system into this one.
		*/
		void CopyPatches(System const& other);


		/**
		\brief Query whether a system is patched.
		*/
		bool IsPatched() const
		{
			return is_patched_;
		}


		template <typename T>
		Vec<T> RescalePointToFitPatch(Vec<T> const& x) const
		{
			return patch_.RescalePoint(x);
		}


		template<typename T>
		void RescalePointToFitPatchInPlace(Vec<T> & x) const
		{
			patch_.RescalePointToFitInPlace(x);
		}
		/**
		 \brief Overloaded operator for printing to an arbirtary out stream.
		 */
		friend std::ostream& operator <<(std::ostream& out, const System & s);


		/**
		 \brief Clear the entire structure of variables in a system.  Reconstructing it is up to you.
		*/
		void ClearVariables();


		/**
		 \brief Copy the entire structure of variables from within one system to another.

		  This copies everything -- ungrouped variables, variable groups, homogenizing variables, the path variable, the ordering of the variables.

		  \param other Another system from which to copy the variable structure.  

		  This operation does NOT affect the functions in any way.  It is up to the user to make the functions actually depend on these variables.
		*/ 
		void CopyVariableStructure(System const& other);
		

		/**
		\brief Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.
		
		\throws std::runtime_error, if the patches are not compatible.  The patches must be either the same, absent, or present in one system.  They propagate to the resulting system.
		*/
		System operator+=(System const& rhs);

		/**
		\brief Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.

		\see The += operator for System also.
		*/
		friend System operator+(System lhs, System const& rhs);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		System operator*=(Nd const& N);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend System operator*(System s, Nd const&  N);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend System operator*(Nd const&  N, System const& s);
	private:

		/**
		\brief Get the sizes according to the FIFO ordering.
		*/
		std::vector<unsigned> VariableGroupSizesFIFO() const;

		/**
		\brief Get the sizes according to the Affine-Homogeneous-Ungrouped ordering.
		*/
		std::vector<unsigned> VariableGroupSizesAffHomUng() const;

	

		/**
		\brief Set up patches automatically for a system using the FIFO ordering.
		*/
		void AutoPatchFIFO();


		/**
		\brief Set up patches automatically for a system using the Affine-Homogeneous-Ungrouped ordering.
		*/
		void AutoPatchAffHomUng();


		/**
		\brief Dehomogenize a point according to the FIFO variable ordering.
	
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

		\see FIFOVariableOrdering
		*/
		template<typename T>
	    Vec<T> DehomogenizePointFIFOOrdering(Vec<T> const& x) const
        {
        	#ifndef BERTINI_DISABLE_ASSERTS
        	assert(ordering_==OrderingChoice::FIFO && "calling FIFO dehomogenize, but system is set to use a different ordering.");
        	assert(homogenizing_variables_.size()==0 || homogenizing_variables_.size()==NumVariableGroups() && "must have either 0 homogenizing variables, or the number of homogenizing variables must match the number of affine variable groups.");
        	#endif

        	bool is_homogenized = homogenizing_variables_.size()!=0;
        	Vec<T> x_dehomogenized(NumNaturalVariables());

        	unsigned affine_group_counter = 0;
        	unsigned hom_group_counter = 0;
        	unsigned ungrouped_variable_counter = 0;

        	unsigned hom_index = 0; // index into x, the point we are dehomogenizing
        	unsigned dehom_index = 0; // index into x_dehomogenized, the point we are computing

    		for (auto iter : time_order_of_variable_groups_)
    		{
    			switch (iter){
    				case VariableGroupType::Affine:
    				{
    					if (is_homogenized)
    					{
	    					auto h = x(hom_index++);
	    					for (unsigned ii = 0; ii < variable_groups_[affine_group_counter].size(); ++ii)
	    						x_dehomogenized(dehom_index++) = x(hom_index++) / h;
	    					affine_group_counter++;
	    				}
	    				else
	    				{
	    					for (unsigned ii = 0; ii < variable_groups_[affine_group_counter].size(); ++ii)
	    						x_dehomogenized(dehom_index++) = x(hom_index++);
	    				}
    					break;
    				}
    				case VariableGroupType::Homogeneous:
    				{
    					for (unsigned ii = 0; ii < hom_variable_groups_[hom_group_counter].size(); ++ii)
    						x_dehomogenized(dehom_index++) = x(hom_index++);
    					break;
    				}
    				case VariableGroupType::Ungrouped:
    				{
    					x_dehomogenized(dehom_index++) = x(hom_index++);
    					ungrouped_variable_counter++;
    					break;
    				}
    				default:
    				{
    					throw std::runtime_error("unacceptable VariableGroupType in FIFOVariableOrdering");
    				}
    			}
    		}

    		return x_dehomogenized;
        }



	    /**
		\brief Dehomogenize a point according to the AffHomUng variable ordering.
		
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

		\see AffHomUngVariableOrdering
		*/
	    template<typename T>
	    Vec<T> DehomogenizePointAffHomUngOrdering(Vec<T> const& x) const
        {
        	#ifndef BERTINI_DISABLE_ASSERTS
        	assert(ordering_==OrderingChoice::AffHomUng && "calling Dehomogenize using Affine-Homogeneous-Ungrouped ordering, but system is in FIFO mode.");
        	#endif

        	Vec<T> x_dehomogenized(NumNaturalVariables());

        	unsigned hom_index = 0; // index into x, the point we are dehomogenizing
        	unsigned dehom_index = 0; // index into x_dehomogenized, the point we are computing

        	for (auto iter=variable_groups_.begin(); iter!=variable_groups_.end(); iter++)
    		{
    			auto h = x(hom_index++);
    			for (unsigned ii = 0; ii < iter->size(); ++ii)
    				x_dehomogenized(dehom_index++) = x(hom_index++) / h;
    		}

    		for (auto iter=hom_variable_groups_.begin(); iter!=hom_variable_groups_.end(); iter++)
    			for (unsigned ii = 0; ii < iter->size(); ++ii)
    				x_dehomogenized(dehom_index++) = x(hom_index++);
    		
    		for (unsigned ii = 0; ii < ungrouped_variables_.size(); ++ii)
    			x_dehomogenized(dehom_index++) = x(hom_index++);

    		return x_dehomogenized;
        }

	    /**
		 Puts together the ordering of variables, and stores it internally.
	    */
	    void ConstructOrdering() const;


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
		
		class Patch patch_; ///< Patch on the variable groups.  Assumed to be in the same order as the time_order_of_variable_groups_ if the system uses FIFO ordering, or in same order as the AffHomUng variable groups if that is set.
		bool is_patched_;	///< Indicator of whether the system has been patched.

		mutable std::vector< Jac > jacobian_; ///< The generated functions from differentiation.  Created when first call for a Jacobian matrix evaluation.
		mutable bool is_differentiated_; ///< indicator for whether the jacobian tree has been populated.


		std::vector< VariableGroupType > time_order_of_variable_groups_;

		mutable std::tuple< Vec<dbl>, Vec<mpfr> > current_variable_values_;

		mutable VariableGroup variable_ordering_; ///< The assembled ordering of the variables in the system.
		mutable bool have_ordering_;
		OrderingChoice ordering_;

		mutable unsigned precision_; ///< the current working precision of the system 


		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {

			ar & ungrouped_variables_;
			ar & variable_groups_;
			ar & hom_variable_groups_;
			ar & homogenizing_variables_;

			ar & time_order_of_variable_groups_;

			ar & have_path_variable_;
			ar & path_variable_;			

			ar & ordering_;
			ar & have_ordering_;
			ar & variable_ordering_;

			ar & implicit_parameters_;
			ar & explicit_parameters_;

			ar & constant_subfunctions_;
			ar & subfunctions_;
			ar & functions_;

			ar & is_differentiated_;
			ar & jacobian_;

			ar & precision_;
			ar & is_patched_;
			ar & patch_;
		}

	};




	




}









#endif // for the ifndef include guards



