//This file is part of Bertini 2.
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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file bertini2/system/system.hpp 

\brief Provides the bertini::System class.
*/

#ifndef BERTINI_SYSTEM_HPP
#define BERTINI_SYSTEM_HPP

#include <assert.h>
#include <vector>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/std_variant.hpp>
#include <boost/type_index.hpp>

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"


#include "bertini2/function_tree.hpp"
#include "bertini2/system/patch.hpp"

#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/blocks/block.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

namespace bertini {

	// (included above) so the evaluation blocks can see them.

	class Slice;  // system/slice.hpp -- a thin wrapper over a LinearFormsBlock; System::Slices()
	              // recovers them from a slice-derived system without binding the block variant.

	/**
	\brief The fundamental polynomial system class for Bertini2.
	
	 The fundamental polynomial system class for Bertini2.

	 Other System types are derived from this, but this class is not abstract.
	 */
	class System{
	
	public:
		// a few local using statements to reduce typing etc.
		using NE = std::shared_ptr<node::NamedExpression>;
		using Var = std::shared_ptr<node::Variable>;
		using Nd = std::shared_ptr<node::Node>;

		/**
		\brief The default constructor for a system.
		*/
		System() : have_path_variable_(false), is_patched_(false), is_differentiated_(false), have_ordering_(false), precision_(DefaultPrecision())
		{}

		/**
		\brief Construct a system from a list of functions.

		The functions are added to the system, the variables appearing in them are
		automatically discovered (see node::GatherVariables), and those variables are
		placed into a single affine variable group (ordered alphabetically by name).
		This is a convenience for programmatically building a system without having to
		assemble the variable group by hand.

		\param functions The functions which define the system.
		*/
		explicit
		System(std::vector<Nd> const& functions);

		/**
		\brief The copy operator, creates a system from a string using the Bertini parser for Bertini classic syntax.
		*/
		explicit
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
		System& operator=(const System & other);

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

		\param new_precision The new precision, in digits, to work in.  This only affects the complex_mp types, not double.  To use low-precision (doubles), use that number type in the templated functions.
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
		 \brief Evaluate the system using the previously set variable (and time) values, in place.

		It is up to YOU to ensure that the system's variables (and path variable) has been set prior to this function call.

		\param function_values The vector to write the function values into.
		*/
		template<typename T>
		void EvalInPlace(Vec<T> & function_values) const
		{
			
			if (function_values.size() != static_cast<Eigen::Index>(NumTotalFunctions())) 
			{
				std::stringstream ss;
				ss << "trying to evaluate system in-place, but number length of vector into which to write the values (" << function_values.size() << ") doesn't match number of system user-defined functions plus patches ( " << NumNaturalFunctions() << "+" << NumPatches() << ") = " << NumTotalFunctions() << ").  Use System.NumTotalFunctions() to make the container for in-place evaluation";
				throw std::runtime_error(ss.str());
			}

			if (!is_differentiated_)
				Differentiate();   // syncs + (for the polynomial block) compiles the SLP

			EvalBlocksInPlace<T>(function_values);
			if (IsPatched())
				patch_.EvalInPlace(function_values,
				                   std::get<Vec<T> >(current_variable_values_));
			CoerceBlockOutputPrecision(function_values);
		}
		
		
		
		
		
		/**
		 \brief Evaluate the system using the previously set variable (and time) values, creating vector of function values.
		 
		 It is up to YOU to ensure that the system's variables (and path variable) has been set prior to this function call.
		 
		 \return The function values of the system
		 */
		template<typename T>
		Vec<T> Eval() const
		{
			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.
			EvalInPlace(function_values);

			return function_values;
		}


		

		/**
		 \brief Evaluate the system, provided the system has no path variable defined, in place.
		 
		 Causes the current variable values to be set in the system.  Resets the function tree's stored numbers.
		 
		 
		 \throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		 \param function_values The vector to write the function values into.
		 \param variable_values The values of the variables, for the evaluation.
		 */
		template<typename T, typename Derived>
		void EvalInPlace(Vec<T>& function_values, const Eigen::MatrixBase<Derived>& variable_values) const
		{
			static_assert(std::is_same<typename Derived::Scalar,T>::value,"scalar types must match");

			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
			{
				std::stringstream ss;
				ss << "trying to evaluate system, but number of input variables (" << variable_values.size() << ") doesn't match number of system variables (" << NumVariables() << ").";
				throw std::runtime_error(ss.str());
			}
			if (have_path_variable_)
				throw std::runtime_error("not using a time value for evaluation of system, but path variable IS defined.");
			
			SetVariables(variable_values.eval());
			EvalInPlace(function_values);
		}
		
		
		
		
		/**
		\brief Evaluate the system, provided the system has no path variable defined.

		Causes the current variable values to be set in the system.  Resets the function tree's stored numbers.  


		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename Derived>
		typename Derived::PlainObject Eval(const Eigen::MatrixBase<Derived>& variable_values) const
		{
			typedef typename Derived::Scalar T;

			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.
			EvalInPlace(function_values, variable_values);
			return function_values;

		}

		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values) const
		{
			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.
			EvalInPlace(function_values, variable_values);
			return function_values;
		}

		
		

		
		/**
		 Evaluate the system, provided a path variable is defined for the system, in place.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		 
		 \param function_values The vector to write the function values into.
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \todo The Eval() function for systems has the unfortunate side effect of resetting constant functions.  Modify the System class so that only certain parts of the tree get reset.
		 */
		template<typename Derived, typename T>
		void EvalInPlace(Vec<T> & function_values, const Eigen::MatrixBase<Derived>& variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");

			SetVariables(variable_values.eval());
			SetPathVariable(path_variable_value);


			EvalInPlace(function_values);
		}

		
		
		
		




		/**
		 Evaluate the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \todo The Eval() function for systems has the unfortunate side effect of resetting constant functions.  Modify the System class so that only certain parts of the tree get reset.
		 */
		template<typename Derived, typename T>
		Vec<T> Eval(const Eigen::MatrixBase<Derived>& variable_values, const T & path_variable_value) const
		{
			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.
			EvalInPlace(function_values, variable_values, path_variable_value);
			return function_values;
		}


		template<typename T>
		Vec<T> Eval(const Vec<T> & variable_values, const T & path_variable_value) const
		{
			Vec<T> function_values(NumTotalFunctions()); // create vector with correct number of entries.
			EvalInPlace(function_values, variable_values, path_variable_value);
			return function_values;
		}
		
		
		
		/**
		 \brief Evaluate the Jacobian matrix of the system, using the previous space and time values, in place.

		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		This is analogous to J = sys.JacobianInPlace();

		The input matrix must have the correct size already
		*/
		template <typename T>
		void JacobianInPlace(Mat<T> & J) const
		{
		

			if(J.rows() != static_cast<Eigen::Index>(NumTotalFunctions()) || J.cols() != static_cast<Eigen::Index>(NumVariables()))
			{
				throw std::runtime_error("trying to evaluate jacobian of system in place, but input J doesn't have right number of columns or rows");
			}
			
			if (!is_differentiated_)
				Differentiate();

			JacobianBlocksInPlace<T>(J);
			if (IsPatched())
				patch_.JacobianInPlace(J, std::get<Vec<T> >(current_variable_values_));
			CoerceBlockOutputPrecision(J);
		}

		
		
		
		
		


		/**
		Evaluate the Jacobian matrix of the system, using the previous space and time values.

		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		*/
		template<typename T>
		Mat<T> Jacobian() const
		{

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J);

			return J;
		}

		

		
		/**
		 Evaluate the Jacobian matrix of the system, provided the system has no path variable defined.
		 
		 \throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		 \param J The matrix to write the Jacobian into.
		 \param variable_values The values of the variables, for the evaluation.
		 */
		template<typename T>
		void JacobianInPlace(Mat<T> & J, const Vec<T> &  variable_values) const
		{

			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");
			
			if (HavePathVariable())
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");
			
			SetAndReset(variable_values);
			
			JacobianInPlace(J);
		}

		

		
		
		
		/**
		Evaluate the Jacobian matrix of the system, provided the system has no path variable defined.

		\throws std::runtime_error, if a path variable IS defined, but you didn't pass it a value.  Also throws if the number of variables doesn't match.
		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values) const
		{
			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (HavePathVariable())
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J,variable_values);
			return J;
		}


		
		
		/**
		 Evaluate the Jacobian of the system, provided a path variable is defined for the system, in place.
		 
		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.

		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		 \param J The matrix to write the Jacobian into.
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.
		 */
		template<typename Derived, typename T>
		void JacobianInPlace(Mat<T> & J, const Eigen::MatrixBase<Derived> & variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");
			
			if (!HavePathVariable())
				throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");
			
			SetVariables(variable_values.eval());
			SetPathVariable(path_variable_value);
			JacobianInPlace(J);
		}

		
		


		/**
		 Evaluate the Jacobian of the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \return The Jacobian matrix.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		 */
		template<typename Derived, typename T>
		Mat<T> Jacobian(const Eigen::MatrixBase<Derived> & variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (!HavePathVariable())
				throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J,variable_values, path_variable_value);
			return J;
		}


		template<typename T>
		Mat<T> Jacobian(const Vec<T> & variable_values, const T & path_variable_value) const
		{
			if (variable_values.size()!=static_cast<Eigen::Index>(NumVariables()))
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (!HavePathVariable())
				throw std::runtime_error("using a time value for computation of jacobian, but no path variable is defined.");

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J,variable_values,path_variable_value);
			return J;
		}

		
		/**
		\brief Compute the time-derivative of a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename Derived, typename T>
		void TimeDerivativeInPlace(Vec<T> & ds_dt, 
							const Eigen::MatrixBase<Derived> & variable_values, 
							const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			SetVariables(variable_values.eval());
			SetPathVariable(path_variable_value);
			TimeDerivativeInPlace(ds_dt);
		}

		
		
		
		
		
		/**
		\brief Compute the time-derivative of a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename Derived, typename T>
		Vec<T> TimeDerivative(const Eigen::MatrixBase<Derived> & variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			Vec<T> ds_dt(NumTotalFunctions());
			TimeDerivativeInPlace(ds_dt, variable_values, path_variable_value);
			return ds_dt;
		}
		







		template<typename Derived, typename T>
		void TimeDerivativeInPlace(Vec<T> & ds_dt, 
							const Eigen::MatrixBase<Derived> & variable_values) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			SetVariables(variable_values.eval());
			TimeDerivativeInPlace(ds_dt);
		}

		
		
		
		
		
		/**
		\brief Compute the time-derivative of a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename Derived, typename T>
		Vec<T> TimeDerivative(const Eigen::MatrixBase<Derived> & variable_values) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			Vec<T> ds_dt(NumTotalFunctions());
			TimeDerivativeInPlace(ds_dt, variable_values);
			return ds_dt;
		}


		template<typename T>
		void TimeDerivativeInPlace(Vec<T> & ds_dt) const
		{

			if(ds_dt.size() < static_cast<Eigen::Index>(NumNaturalFunctions()))
			{
				std::stringstream ss;
				ss << "trying to evaluate system in place, but number of input functions (" << ds_dt.size() << ") doesn't match number of system functions (" << NumNaturalFunctions() << ").";
				throw std::runtime_error(ss.str());
			}

			if (!HavePathVariable())
				throw std::runtime_error("computing time derivative of system with no path variable defined");

			if (!is_differentiated_)
				Differentiate();

			TimeDerivBlocksInPlace<T>(ds_dt);
			// the patch doesn't move with time.  derivatives 0.
			if (IsPatched())
				for (size_t ii = 0; ii < NumTotalVariableGroups(); ++ii)
					ds_dt(static_cast<Eigen::Index>(ii + NumNaturalFunctions())) = T(0);
			CoerceBlockOutputPrecision(ds_dt);
		}

		
		
		
		
		
		/**
		\brief Compute the time-derivative of a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename T>
		Vec<T> TimeDerivative() const
		{
			if (!HavePathVariable())
				throw std::runtime_error("computing time derivative of system with no path variable defined");

			Vec<T> ds_dt(NumTotalFunctions());
			TimeDerivativeInPlace(ds_dt);
			return ds_dt;
		}
		
		/**
		Homogenize the system, adding new homogenizing variables for each VariableGroup defined for the system.

		\throws std::runtime_error, if the system is not polynomial, has a mismatch on the number of homogenizing variables and the number of variable groups (this would result from a partially homogenized system), or the homogenizing variable names somehow get screwed up by having duplicates.
		*/
		void Homogenize();

		/**
		Homogenize the system reusing externally-supplied homogenizing variables -- one per affine
		variable group, in group order -- instead of minting fresh ones.  Used by
		`RandomizationBlock` so a wrapped operand system shares the owning system's homogenizing
		variables (the h-power deficit factors and the operand's functions must live in the same
		homogeneous coordinates).  Otherwise behaves exactly like `Homogenize()` on a fresh system.

		\throws std::runtime_error if the system is non-polynomial, already homogenized, or the
		number of supplied homogenizing variables does not equal the number of affine variable groups.
		*/
		void Homogenize(VariableGroup const& provided_hom_vars);

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
		size_t NumNaturalFunctions() const;

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
		 Get the homogenizing variables, one per affine variable group (in group order), or
		 an empty container if the system is not homogenized.
		 */
		VariableGroup const& HomogenizingVariables() const { return homogenizing_variables_; }

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

		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.
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
			if (new_values.size()!= static_cast<Eigen::Index>(NumVariables()))
			{
				throw std::runtime_error("variable vector of different length from system-owned variables in SetVariables");
			}

			#ifndef BERTINI_DISABLE_PRECISION_CHECKS
				// A system with no variables (a constant) has an empty point: there is no
				// precision to read from it, so skip the check.
				if (new_values.size() > 0)
				{
					if constexpr (!std::is_same<T,complex_dbl>::value) {
						if (Precision(new_values) != this->precision())
							throw std::runtime_error("precision of input point in SetVariables (" + std::to_string(Precision(new_values)) + ") must match the precision of the system (" + std::to_string(this->precision()) + ").");
					}
				}
			#endif

			// Blocks are value-in: the polynomial block feeds this stored vector into its SLP, the
			// structured blocks compute on it directly, and the patch reads it too (see
			// EvalBlocksInPlace).  The shared Variable nodes are no longer written during evaluation
			// (ADR-0027), so the node DAG stays read-only across threads.
			std::get<Vec<T> >(current_variable_values_) = new_values;
		}



		/**
		 Set the current value of the path variable.
		
		 \throws std::runtime_error, if a path variable is not defined.
		 \tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		 \param new_value The new updated values for the path variable.
		 */
		template<typename T>
		void SetPathVariable(T const& new_value) const
		{
			if (!have_path_variable_)
				throw std::runtime_error("trying to set the value of the path variable, but one is not defined for this system");

			// Store the path value in the System's own per-thread buffer, NOT the shared node.
			// Blocks are value-in and receive the path value as an argument (see EvalBlocksInPlace /
			// CurrentPathValue), so the path-variable node is never read during evaluation (ADR-0027).
			std::get<T>(current_path_value_) = new_value;
		}


		// Stage the system's current point (variables, and optionally the path value) for a
		// subsequent Eval/Jacobian.  (Formerly also reset the function-tree node caches; node-level
		// evaluation is gone, so there is nothing to reset -- the SLP carries its own state.)
		template<typename T>
		void SetAndReset(Vec<T> const& new_space, T const& new_time) const
		{
			SetVariables(new_space);
			SetPathVariable(new_time);
		}

		template<typename T>
		void SetAndReset(Vec<T> const& new_space) const
		{
			SetVariables(new_space);
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
		 \brief Replace the entire variable-group structure of the system.

		 Clears all existing affine, homogeneous, and ungrouped variables (and any
		 homogenizing variables introduced by a prior Homogenize), then installs the
		 supplied groups as the system's affine variable groups, in order.  The path
		 variable, if any, is preserved.

		 \param groups The affine variable groups to install.
		 */
		void SetVariableGroups(std::vector<VariableGroup> const& groups);


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
		void AddParameter(NE const& F);



		/**
		 Add a function to the system, as a bare expression.

		 \param F The function to add.
		 */
		void AddFunction(Nd const& F);


		/**
		 Add some functions to the system.

		 \param F The functions to add.
		 */
		void AddFunctions(std::vector<Nd> const& F);

		/**
		\brief Append an evaluation block (products-of-linears, blend, ...).

		A block-composed system evaluates its blocks (in the order added) instead of the
		function-tree / SLP functions; blocks contribute the leading "natural" rows, with
		any patch appended after, exactly as for a classic system.
		*/
		void AddBlock(Block b) { blocks_.push_back(std::move(b)); }

		/// Read-only access to the system's evaluation blocks (in evaluation order).  For
		/// introspection / testing -- e.g. reading a RandomizationBlock's matrix; the variant
		/// itself stays out of the Python surface.
		std::vector<Block> const& Blocks() const { return blocks_; }

		/// \brief Recover the linear-form slices embedded in this system, one per LinearFormsBlock
		/// (empty if the system has none).  Lets a user back out the slice structure of a system
		/// built from a slice -- the Python-facing alternative to exposing the Block variant.
		/// Defined in system.cpp (needs the full Slice type).
		std::vector<Slice> Slices() const;

		/// Remove the polynomial block (the System's natural functions), leaving any structured
		/// blocks and the variable structure / patch intact.  Used to turn a copy of a System
		/// into a homotopy shell whose rows come from a blend block rather than its own
		/// functions (FormHomotopy): `h = target; h.ClearFunctions(); h.AddBlock(blend);`.
		void ClearFunctions()
		{
			for (auto it = blocks_.begin(); it != blocks_.end(); )
			{
				if (std::holds_alternative<blocks::PolynomialBlock>(*it))
					it = blocks_.erase(it);
				else
					++it;
			}
			InvalidateDifferentiation();
		}

		/// \brief Whether this system is evaluated from blocks rather than the function tree.
		bool HasBlocks() const { return !blocks_.empty(); }

		/// Does the system have any non-polynomial (structured) block -- products-of-linears,
		/// linear-forms, blend?  Every system has a PolynomialBlock for its functions after the
		/// fold, so HasBlocks() is no longer the right test for "needs whole-System blending";
		/// this is.  (e.g. an MHom start system has a products block; a total-degree start does not.)
		bool HasStructuredBlocks() const
		{
			for (auto const& b : blocks_)
				if (!std::holds_alternative<blocks::PolynomialBlock>(b))
					return true;
			return false;
		}

		/// \brief Remove all evaluation blocks (the system reverts to its function-tree
		/// functions).  Mainly for testing the block path against the function-tree path.
		void ClearBlocks() { blocks_.clear(); }

		/// \brief Build an equivalent **pure function-tree** System: every block's functions
		/// expressed as function-tree nodes, gathered into a single PolynomialBlock, with the
		/// same variables, path variable, and patch (the patch is reused, not re-expressed).
		///
		/// This is a verification / interop oracle — the block path exists for performance and
		/// precision control, so this is NOT a replacement for block evaluation.  It lets the
		/// block-composed evaluation be cross-checked against the function-tree path
		/// (eval / Jacobian must agree).  Scoped to the current block types; a block that cannot
		/// be expanded throws.
		System ExpandToFunctionTree() const;

		/// \brief The system's natural (pre-patch) functions as function-tree expression nodes,
		/// expanding any structured block.  Used by ExpandToFunctionTree and, recursively, by
		/// BlendBlock expansion (a blend is sum_i c_i(t) * operand_i, each operand expanded).
		std::vector<Nd> NaturalFunctionsAsNodes() const;


		/**
		\brief Randomize an overdetermined system down to a square one, returning a NEW system.

		An overdetermined system (N natural functions, n variables, N > n) is replaced by n generic
		combinations whose isolated solutions still contain this system's -- the standard squaring-up
		so the isolated solutions can be found by homotopy continuation (solve the square result,
		then discard the extraneous solutions by re-evaluating this system).

		The returned system carries a single `RandomizationBlock` wrapping a copy of this system; the
		combination is `g_i = sum_j R_ij f_j` (the block applies the homogenizing-variable powers
		needed when the functions differ in degree).  **This system is not mutated** -- any internal
		degree sorting happens on the copy.

		Auto form: for a single affine variable group the functions are sorted by descending degree
		and `R = [I | C]` (C random), giving the optimal total-degree path count (the product of the
		n largest degrees); for several variable groups a dense random `R` is used with a common
		target multidegree.

		\throws std::runtime_error if the system is underdetermined (fewer functions than variables).
		*/
		System Randomize() const;

		/**
		\brief Randomize using a caller-supplied (exact) coefficient matrix R, leaving the functions
		in their current order.  R has one row per desired randomized function and one column per
		natural function of this system.  Returns a NEW system; this one is not mutated.

		\throws std::runtime_error if R's column count does not equal this system's natural-function count.
		*/
		System Randomize(Mat<complex_mp> const& R) const;

		/**
		\brief The randomization matrix R of a system produced by Randomize() (its first
		RandomizationBlock's n x N coefficient matrix).

		\throws std::runtime_error if this system carries no randomization block.
		*/
		Mat<complex_mp> RandomizationMatrix() const;







		/**
		 Add a constant function to the system.  Constants must not depend on anything which can vary -- they're constant!

		 \param C The constant to add.
		 */
		void AddConstant(NE const& C);




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

		const Var& GetPathVariable() const;

		/**
		 Order the variables, by the order in which the groups were added.

		 This function returns the variables in First In First Out (FIFO) ordering.

		 Homogenizing variables precede affine variable groups, so that groups always are grouped together.

		This function freshly constructs the ordering every time. If you want to use the cached value, \see Variables.

		 \throws std::runtime_error, if there is a mismatch between the number of homogenizing variables and the number of variable_groups.  This would happen if a system is homogenized, and then more stuff is added to it.  
		*/
		VariableGroup VariableOrdering() const;


		/**
		 Get the variables in the problem.

		 Returns the variable ordering, constructing it if necessary.
		*/
		const VariableGroup& Variables() const;

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
			return VariableGroupSizesFIFO();
		}

		/**
		\brief Dehomogenize a point, using the variable grouping / structure of the system.
		
		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		\throws std::runtime_error, if there is a mismatch between the number of variables in the input point, and the total number of var
		*/
		template<typename T>
		Vec<T> DehomogenizePoint(Vec<T> const& x) const
			{

				if (x.size()!=static_cast<Eigen::Index>(NumVariables())){
					std::stringstream message;
					message << "dehomogenizing point with incorrect number of coordinates. input has ";
					message << x.size();
					message << " but system expects ";
					message << NumVariables();
					throw std::runtime_error(message.str());
				}

				if (!have_ordering_)
					ConstructOrdering();

				return DehomogenizePointFIFO(x);
			}


		/**
		\brief The infinity norm of a point after dehomogenization.

		This is the single canonical "how big is this point, in user coordinates" measurement.
		An endpoint going to infinity has its dehomogenized coordinates blow up, so this is what
		the endgames test against `Security::max_norm` to detect divergence, and what the
		zero-dim solver tests against `endpoint_finite_threshold` to classify finite/infinite
		endpoints.  Routing both through here keeps those decisions consistent: never compare the
		raw internal (homogenized, on-patch) coordinates, which carry the homogenizing variable
		and patch scaling.

		\tparam T the number-type of the point.  Returns the associated real magnitude type.
		*/
		template<typename T>
		auto InfinityNormOfDehomogenized(Vec<T> const& x) const
		{
			return DehomogenizePoint(x).template lpNorm<Eigen::Infinity>();
		}


		/**
		\brief Take a point in user (dehomogenized) coordinates into this system's internal coordinates.

		Two steps: (1) insert the homogenizing coordinate, with value 1, for each affine
		variable group, per the FIFO variable ordering; (2) if the system is patched,
		rescale the result onto the patch.  The result is the projectively-identical
		point expressed in the coordinates the solver works in -- suitable for
		comparison with internal-coordinate solutions, or as a start point for further
		tracking re-using this system's patch.

		This is the inverse of DehomogenizePoint: DehomogenizePoint(HomogenizePoint(p)) == p.
		On an unhomogenized, unpatched system this is the identity.

		\tparam T the number-type.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		\throws std::runtime_error, if there is a mismatch between the number of variables in the input point, and the number of natural variables of the system.
		*/
		template<typename T>
		Vec<T> HomogenizePoint(Vec<T> const& x) const
			{

				if (x.size()!=static_cast<Eigen::Index>(NumNaturalVariables())){
					std::stringstream message;
					message << "homogenizing point with incorrect number of coordinates. input has ";
					message << x.size();
					message << " but system expects ";
					message << NumNaturalVariables();
					throw std::runtime_error(message.str());
				}

				if (!have_ordering_)
					ConstructOrdering();

				auto x_homogenized = HomogenizePointFIFO(x);

				if (IsPatched())
					RescalePointToFitPatchInPlace(x_homogenized);

				return x_homogenized;
			}




		/**
		 \brief Get a function by its index, as a function-tree node.

		 Works for every block type: structured blocks (linear forms, products of linears,
		 randomization, blend) are expanded to their node form on demand, so this is in sync with
		 NumNaturalFunctions().  Throws std::out_of_range if the index is past the end (rather than
		 dereferencing past it -- the old un-checked version segfaulted on slice-derived systems
		 that have no PolynomialBlock).
		*/
		auto Function(unsigned index) const
		{
			auto fns = NaturalFunctionsAsNodes();
			if (index >= fns.size())
				throw std::out_of_range("System::Function index out of range");
			return fns[index];
		}


		/**
		 \brief Get the functions, as function-tree nodes.

		 Expands any structured block (so it agrees with NumNaturalFunctions() and Function(i),
		 even for systems built purely from a slice / start-system block).
		*/
		std::vector<Nd> GetNaturalFunctions() const
		{
			return NaturalFunctionsAsNodes();
		}



		/**
		 Get the affine variable groups in the problem.
		*/
		auto VariableGroups() const
		{
			return variable_groups_;
		}

		/**
		 Get the homogeneous (projective) variable groups in the problem.
		*/
		auto HomVariableGroups() const
		{
			return hom_variable_groups_;
		}


		const Var& PathVariable() const
		{
			return path_variable_;
		}


		/**
		 \brief Get the space derivatives

		 These are as computed using the Derivatives method, not Jacobian nodes.  They are stored in column-major order, so stride over variables first.

		 ```
		 for (int jj = 0; jj < num_vars; ++jj)
			for (int ii = 0; ii < num_functions; ++ii)
				space_derivatives_[ii+jj*num_functions] = functions_[ii]->Differentiate(vars[jj]);
		 ```
		 */
		std::vector< Nd > GetSpaceDerivatives() const;

		/**
		 \brief Get the time derivatives of all functions

		 These are as computed using the Derivatives method, not Jacobian nodes.  Stored in order of functions.
		 */
		std::vector< Nd > GetTimeDerivatives() const;


		//////////////////////
		//
		//  Functions involving coefficients of the system
		//
		///////////////////////


		/**
		Compute an estimate of an upper bound of the absolute values of the coefficients in the system.
		
		\param num_evaluations The number of times to compute this estimate.  Default is 1.
		\returns An upper bound on the absolute values of the coefficients.
		*/
		template <typename NumT>
		typename Eigen::NumTraits<NumT>::Real CoefficientBound(unsigned num_evaluations = 1) const;
		


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
			AutoPatchFIFO();
		}

		/**
		\brief Copy the patches from another system into this one.
		*/
		void CopyPatches(System const& other);


		Patch GetPatch() const
		{
			return patch_;
		}

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
		\brief Human-facing description of the system, block by block.

		Each block describes the rows it owns (labelled `f_k` in function-vector order).  `verbose ==
		false` (the default, used by `operator<<` / Python `str`) shows placeholder symbols for the
		structured blocks' matrices/coefficients; `verbose == true` reveals the actual numbers and the
		underlying functions of randomization / blend blocks.  This is for reading, not re-parsing.
		*/
		void Describe(std::ostream& out, bool verbose = false) const;

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
		
		// The Please/Dont AssumeUniformPrecision family was removed: the setter had
		// ignored its argument (always storing false) for ages, so the early-out in
		// System::precision() it was meant to enable was dead code, and skipping the
		// propagation is unsound anyway (e.g. the SLP can be at a different precision
		// than precision_ claims).  precision() now always propagates.

		/**
		\brief Simplify the functions contained in the system.

		\note This may change any nodes on which the system depends.
		*/
		void SimplifyFunctions();

		/**
		\brief Simplify the derivatives / jacobian / etc contained in the system.

		\note This may change any nodes on which the system depends.
		*/
		void SimplifyDerivatives() const;

		/**
		\brief Simplify as many aspects of the system as possible.  

		\note This may change any nodes on which the system depends.
		*/
		void Simplify();





		/**
		\brief Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.
		
		\throws std::runtime_error, if the patches are not compatible.  The patches must be either the same, absent, or present in one system.  They propagate to the resulting system.
		*/
		System& operator+=(System const& rhs);

		/**
		\brief Add two systems together.

		\throws std::runtime_error, if the systems are not of compatible size -- either in number of functions, or variables.  Does not check the structure of the variables, just the numbers.

		\see The += operator for System also.
		*/
		friend const System operator+(System lhs, System const& rhs);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		System& operator*=(Nd const& N);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend const System operator*(System s, Nd const&  N);

		/**
		\brief Multiply a system by an arbitrary node.  

		Can be used for defining a coupling of a target and start system through a path variable.  Does not affect path variable declaration, or anything else.  It is up to you to ensure the system depends on this node properly.
		*/
		friend const System operator*(Nd const&  N, System const& s);
	private:

		/// Shared back end of the Randomize overloads: given the overdetermined operand (already a
		/// copy, sorted or not) and a finished coefficient matrix, compute the per-row target
		/// multidegrees, build the RandomizationBlock, and return a new system carrying it.
		System AssembleRandomized(std::shared_ptr<System> operand, Mat<complex_mp> coefficients) const;

		/**
		\brief Get the sizes according to the FIFO ordering.
		*/
		std::vector<unsigned> VariableGroupSizesFIFO() const;


		/**
		\brief Set up patches automatically for a system using the FIFO ordering.
		*/
		void AutoPatchFIFO();



		/**
		\brief Dehomogenize a point according to the FIFO variable ordering.
	
		\tparam T the number-type for return.  Probably complex_dbl=std::complex<double>, or complex_mp=bertini::complex_mp.

		\see FIFOVariableOrdering
		*/
		template<typename T>
		Vec<T> DehomogenizePointFIFO(Vec<T> const& x) const
		{
			#ifndef BERTINI_DISABLE_ASSERTS
			assert(homogenizing_variables_.size()==0 || homogenizing_variables_.size()==NumVariableGroups() && "must have either 0 homogenizing variables, or the number of homogenizing variables must match the number of affine variable groups.");
			#endif

			bool is_homogenized = homogenizing_variables_.size()!=0;
			Vec<T> x_dehomogenized(NumNaturalVariables());

			unsigned affine_group_counter = 0;
			unsigned hom_group_counter = 0;

			unsigned hom_index = 0; // index into x, the point we are dehomogenizing
			unsigned dehom_index = 0; // index into x_dehomogenized, the point we are computing

			for (auto& iter : time_order_of_variable_groups_)
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
						hom_group_counter++; // was missing; mattered only for multiple hom groups of differing sizes
						break;
					}
					case VariableGroupType::Ungrouped:
					{
						x_dehomogenized(dehom_index++) = x(hom_index++);
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
		\brief FIFO-ordering implementation of HomogenizePoint's first step: insert
		the homogenizing coordinate, with value 1, at each affine group's slot.
		Homogeneous groups and ungrouped variables pass through.  Patch rescaling is
		the caller's job.
		*/
		template<typename T>
		Vec<T> HomogenizePointFIFO(Vec<T> const& x) const
		{
			#ifndef BERTINI_DISABLE_ASSERTS
			assert(homogenizing_variables_.size()==0 || homogenizing_variables_.size()==NumVariableGroups() && "must have either 0 homogenizing variables, or the number of homogenizing variables must match the number of affine variable groups.");
			#endif

			bool is_homogenized = homogenizing_variables_.size()!=0;
			if (!is_homogenized)
				return x;

			Vec<T> x_homogenized(NumVariables());

			unsigned affine_group_counter = 0;
			unsigned hom_group_counter = 0;

			unsigned dehom_index = 0; // index into x, the user-coordinates point
			unsigned hom_index = 0; // index into x_homogenized, the point we are computing

			for (auto& iter : time_order_of_variable_groups_)
			{
				switch (iter){
					case VariableGroupType::Affine:
					{
						x_homogenized(hom_index++) = T(1);
						for (unsigned ii = 0; ii < variable_groups_[affine_group_counter].size(); ++ii)
							x_homogenized(hom_index++) = x(dehom_index++);
						affine_group_counter++;
						break;
					}
					case VariableGroupType::Homogeneous:
					{
						for (unsigned ii = 0; ii < hom_variable_groups_[hom_group_counter].size(); ++ii)
							x_homogenized(hom_index++) = x(dehom_index++);
						hom_group_counter++;
						break;
					}
					case VariableGroupType::Ungrouped:
					{
						x_homogenized(hom_index++) = x(dehom_index++);
						break;
					}
					default:
					{
						throw std::runtime_error("unacceptable VariableGroupType in HomogenizePointFIFO");
					}
				}
			}

			return x_homogenized;
		}

		// --- the System's polynomial block (its functions live here after the fold) ---

		/// Get the System's PolynomialBlock, creating an (empty) one in blocks_ if none exists.
		blocks::PolynomialBlock& PolyBlock()
		{
			for (auto& b : blocks_)
				if (auto* p = std::get_if<blocks::PolynomialBlock>(&b))
					return *p;
			blocks_.emplace_back(blocks::PolynomialBlock{});
			return std::get<blocks::PolynomialBlock>(blocks_.back());
		}

		/// Find the System's PolynomialBlock, or nullptr if it has none (e.g. a pure
		/// structured system, or a freshly-constructed one with no functions yet).
		blocks::PolynomialBlock const* PolyBlockPtr() const
		{
			for (auto const& b : blocks_)
				if (auto* p = std::get_if<blocks::PolynomialBlock>(&b))
					return p;
			return nullptr;
		}

		/// The polynomial block's functions (an empty list if there is no polynomial block).
		std::vector<Nd> const& PolyFunctions() const
		{
			static const std::vector<Nd> none;
			auto* p = PolyBlockPtr();
			return p ? p->Functions() : none;
		}

		/// Mark the blocks as needing (re)differentiation after a structural change.
		void InvalidateDifferentiation() const
		{
			is_differentiated_ = false;
			if (auto* p = PolyBlockPtr())
				p->Invalidate();
		}

		/// Push the System-owned context (variable ordering, path variable, auto-simplify) into
		/// the PolynomialBlock before it differentiates/evaluates.  The block's setters are
		/// idempotent (they only invalidate on real change), so this is safe to call repeatedly.
		void SyncPolyBlock() const
		{
			if (auto* p = PolyBlockPtr())
			{
				p->SetVariableOrdering(Variables());
				if (have_path_variable_) p->SetPathVariable(path_variable_);
				else                     p->ClearPathVariable();
			}
		}

		/**
		 Puts together the ordering of variables, and stores it internally.
		*/
		void ConstructOrdering() const;

		/// \brief The current path-variable value as type T (zero if no path variable).
		template <typename T>
		T CurrentPathValue() const
		{
			if (have_path_variable_)
				return std::get<T>(current_path_value_);
			return T(0);
		}

		/// \brief Force a block-path evaluation result to the system's working precision.
		///
		/// A block evaluates correctly at its working precision, but the *result* container
		/// (function values / Jacobian / time derivative) is allocated by the caller, often at
		/// whatever the ambient DefaultPrecision happens to be, and Eigen's coefficient-wise
		/// assignment into it preserves the destination entry's precision.  So writing a
		/// 20-digit block value into a result entry that was allocated at, say,
		/// MaxPrecisionAllowed leaves a 20-digit value carried at 1000-digit precision.  The
		/// adaptive tracker then propagates that over-precise value as the path point and the
		/// next System::SetVariables throws (point precision != system precision).  Coercing the
		/// whole result to precision_ here makes a block-composed System honor the contract that
		/// its evaluations come out at its working precision, exactly as the SLP path does.
		/// No-op for double (which carries no precision).
		template <typename Derived>
		void CoerceBlockOutputPrecision(Eigen::MatrixBase<Derived>& result) const
		{
			using Scalar = typename Derived::Scalar;
			if constexpr (!std::is_same<Scalar, complex_dbl>::value)
			{
				using bertini::Precision;
				Precision(result, precision_);
			}
		}

		/// \brief Evaluate the blocks' function values into the leading (natural) rows.
		template <typename T>
		void EvalBlocksInPlace(Vec<T>& function_values) const
		{
			const auto& vars = std::get<Vec<T> >(current_variable_values_);
			const T t = CurrentPathValue<T>();
			Eigen::Index row = 0;
			for (auto const& blk : blocks_)
				std::visit([&](auto const& b){
					const Eigen::Index n = static_cast<Eigen::Index>(b.NumFunctions());
					b.template EvalInPlace<T>(function_values.segment(row, n), vars, t);
					row += n;
				}, blk);
		}

		/// \brief Evaluate the blocks' Jacobian into the leading rows of J.
		template <typename T>
		void JacobianBlocksInPlace(Mat<T>& J) const
		{
			const auto& vars = std::get<Vec<T> >(current_variable_values_);
			const T t = CurrentPathValue<T>();
			Eigen::Index row = 0;
			for (auto const& blk : blocks_)
				std::visit([&](auto const& b){
					const Eigen::Index n = static_cast<Eigen::Index>(b.NumFunctions());
					Mat<T> jb(n, J.cols());                 // blocks write into a contiguous target
					b.template JacobianInPlace<T>(jb, vars, t);
					J.block(row, 0, n, J.cols()) = jb;
					row += n;
				}, blk);
		}

		/// \brief Evaluate the blocks' time-derivative into the leading rows.
		template <typename T>
		void TimeDerivBlocksInPlace(Vec<T>& ds_dt) const
		{
			const auto& vars = std::get<Vec<T> >(current_variable_values_);
			const T t = CurrentPathValue<T>();
			Eigen::Index row = 0;
			for (auto const& blk : blocks_)
				std::visit([&](auto const& b){
					const Eigen::Index n = static_cast<Eigen::Index>(b.NumFunctions());
					b.template TimeDerivInPlace<T>(ds_dt.segment(row, n), vars, t);
					row += n;
				}, blk);
		}


		VariableGroup ungrouped_variables_; ///< ungrouped variable nodes.  Not in an affine variable group, not in a projective group.  Just hanging out, being a variable.
		std::vector< VariableGroup > variable_groups_; ///< Affine variable groups.  When system is homogenized, will have a corresponding homogenizing variable.
		std::vector< VariableGroup > hom_variable_groups_; ///< Homogeneous or projective variable groups.  System SHOULD be homogeneous with respect to these.  

		VariableGroup homogenizing_variables_; ///< homogenizing variables for the variable_groups.  


		bool have_path_variable_ = false; ///< Whether we have the variable or not.
		Var path_variable_; ///< the single path variable for this system.  Sometimes called time.
		
		VariableGroup implicit_parameters_; ///< Implicit parameters.  These don't depend on anything, and will be moved from one parameter point to another by the tracker.  They should be algebraically constrained by some equations.
		std::vector< NE > explicit_parameters_; ///< Explicit parameters.  These should be functions of the path variable only, NOT of other variables.

		// The polynomial path -- functions_, subfunctions_, constant_subfunctions_, their
		// derivatives, the SLP, and eval_method_ -- has been folded into a
		// blocks::PolynomialBlock held in blocks_ (see PolyBlock()/PolyBlockPtr()).  The System
		// is now a thin orchestrator over blocks + variable groups + patch.

		class Patch patch_; ///< Patch on the variable groups.  Assumed to be in the same order as the time_order_of_variable_groups_ if the system uses FIFO ordering, or in same order as the AffHomUng variable groups if that is set.
		bool is_patched_ = false;	///< Indicator of whether the system has been patched.

		mutable bool is_differentiated_ = false; ///< orchestrator flag: have the blocks been differentiated + synced since the last structural change.

		std::vector<Block> blocks_; ///< Evaluation blocks.  Every System has one (a PolynomialBlock for its functions); structured systems (MHom, linear forms) add more.

		std::vector< VariableGroupType > time_order_of_variable_groups_;

		mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > current_variable_values_;
		mutable std::tuple< complex_dbl, complex_mp > current_path_value_{}; ///< per-thread path value (node-free; read by CurrentPathValue, written by SetPathVariable)

		mutable VariableGroup variable_ordering_; ///< The assembled ordering of the variables in the system.
		mutable bool have_ordering_ = false;

		mutable unsigned precision_; ///< the current working precision of the system







		friend class boost::serialization::access;

		/*definition of serialize function*/
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /*version*/){
			ar & ungrouped_variables_;
			ar & variable_groups_;
			ar & hom_variable_groups_;

			ar & homogenizing_variables_;

			ar & have_path_variable_;
			ar & path_variable_;			

			ar & implicit_parameters_;
			ar & explicit_parameters_;

			ar & patch_;
			ar & is_patched_;

			// The polynomial path (functions / subfunctions / derivatives / SLP / eval+deriv
			// methods) now lives inside the PolynomialBlock, which is archived as part of blocks_.
			ar & blocks_;


			// now for the cached / mutable things
			ar & precision_;

			ar & is_differentiated_;

			ar & time_order_of_variable_groups_;

			// if (Archive::is_loading::value == true){
				// have_ordering_ = false;}
			// else
			// {
				ar & have_ordering_;
				ar & variable_ordering_;
			// }


			// if (Archive::is_loading::value == true){
        	// 	std::get<Vec<complex_dbl>>(current_variable_values_).resize(NumVariables());
        	// 	std::get<Vec<complex_mp>>(current_variable_values_).resize(NumVariables());
			// }
			// // next two lines no matter what
			// ar & std::get<Vec<complex_dbl>>(current_variable_values_);
			// ar & std::get<Vec<complex_mp>>(current_variable_values_);

		}


	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};


	/**
	\brief Concatenate two compatible systems.

	Two systems are compatible for concatenation if they have the same variable structure, and if they have the same patch (if patched).

	\param sys1 The top system.
	\param sys2 The bottom system.

	If both patched both must have same patch.  If not both are patched, then the patch will propagate to the returned system. 

	If the two patches have differing variable orderings, the call to Concatenate will throw.
	*/
	System Concatenate(System sys1, System const& sys2);
	


	/**
	\brief Form the gamma-trick straight-line homotopy H = (1-t)*target + gamma*t*start.

	The path variable `t` (named `path_variable_name`) is added to the returned homotopy, tracked
	from t=1 (where H is gamma*start, so its roots are start's solutions) down to t=0 (where H is
	target).  When `start` carries a structured evaluation block (e.g. a products-of-linears start
	system) it cannot be fused by node arithmetic, so the two systems are combined with a
	BlendBlock that evaluates whole Systems; otherwise the node-arithmetic combination is used.
	This is the same construction the zero-dim solver's CloneGiven policy uses internally; it is
	exposed so a user-authored start system can be turned into a trackable homotopy for the
	user-homotopy solve path.

	\param target The target system (H at t=0).
	\param start The start system whose solutions seed the homotopy (H at t=1).
	\param path_variable_name The name to give the homotopy's path variable.
	\param gamma The gamma coefficient (a node).  If null, a random rational gamma is generated.
	*/
	System MakeHomotopy(System const& target, System const& start,
	                    std::string const& path_variable_name = "t",
	                    std::shared_ptr<node::Node> const& gamma = nullptr);

	/**
	\brief Form a homotopy that moves ONLY some rows, leaving the rest fixed and evaluated once.

	The "regeneration" homotopy: `fixed` holds the equations that do not move (the polynomial system
	and any static linear slices); `start_moving` and `end_moving` are systems holding just the rows
	that move, agreeing in function count and variable structure.  The result is

	    H = [ fixed's blocks (unchanged) ;  (1-t)*end_moving + gamma*t*start_moving ]

	with the moving rows produced by a single `BlendBlock` over `end_moving`/`start_moving` and the
	fixed rows kept as their own sibling blocks.  Because the fixed blocks are autonomous, they are
	evaluated exactly once per point and contribute zero to `dH/dt` -- the fixed system is never
	re-evaluated or scaled by the path coefficient as the moving rows slide.  At t=1 the moving rows
	are `gamma*start_moving` (so the start points are roots of `fixed` together with `start_moving`),
	at t=0 they are `end_moving`.

	The fixed rows come first, then the moving rows; build the matching target (for the user-homotopy
	pipeline) as `fixed` concatenated with `end_moving`, and the start points as the roots of `fixed`
	together with `start_moving`.

	\param fixed The equations that do not move (evaluated once per point).
	\param start_moving The moving rows at t=1 (their roots, with `fixed`, are the start points).
	\param end_moving The moving rows at t=0 (the target rows).
	\param path_variable_name The name to give the homotopy's path variable.
	\param gamma The gamma coefficient (a node).  If null, a random rational gamma is generated.
	*/
	System MakeMovingHomotopy(System const& fixed, System const& start_moving, System const& end_moving,
	                          std::string const& path_variable_name = "t",
	                          std::shared_ptr<node::Node> const& gamma = nullptr);

	/**
	\brief Do a deep clone of the system.  This includes the entire structure, variables, etc.  everything.
	*/
	System Clone(System const& sys);

	/**
	\brief Free form function for simplifying systems.
	*/
	void Simplify(System & sys);

	/// \cond INTERNAL
	// Explicit instantiation declarations for the two concrete numeric types.
	// Definitions live in core/src/system/system.cpp.
	// Suppresses re-instantiation of the heavy Eval/Jacobian/Set template bodies
	// (with their eval_method_ switch trees) in every including TU.

	extern template void System::EvalInPlace<complex_dbl>(Vec<complex_dbl>&) const;
	extern template void System::EvalInPlace<complex_mp>(Vec<complex_mp>&) const;

	extern template Vec<complex_dbl> System::Eval<complex_dbl>() const;
	extern template Vec<complex_mp> System::Eval<complex_mp>() const;

	extern template void System::JacobianInPlace<complex_dbl>(Mat<complex_dbl>&) const;
	extern template void System::JacobianInPlace<complex_mp>(Mat<complex_mp>&) const;

	extern template Mat<complex_dbl> System::Jacobian<complex_dbl>() const;
	extern template Mat<complex_mp> System::Jacobian<complex_mp>() const;

	extern template Mat<complex_dbl> System::Jacobian<complex_dbl>(const Vec<complex_dbl>&) const;
	extern template Mat<complex_mp> System::Jacobian<complex_mp>(const Vec<complex_mp>&) const;

	extern template void System::JacobianInPlace<complex_dbl>(Mat<complex_dbl>&, const Vec<complex_dbl>&) const;
	extern template void System::JacobianInPlace<complex_mp>(Mat<complex_mp>&, const Vec<complex_mp>&) const;

	extern template void System::TimeDerivativeInPlace<complex_dbl>(Vec<complex_dbl>&) const;
	extern template void System::TimeDerivativeInPlace<complex_mp>(Vec<complex_mp>&) const;

	extern template Vec<complex_dbl> System::TimeDerivative<complex_dbl>() const;
	extern template Vec<complex_mp> System::TimeDerivative<complex_mp>() const;

	extern template void System::SetVariables<complex_dbl>(const Vec<complex_dbl>&) const;
	extern template void System::SetVariables<complex_mp>(const Vec<complex_mp>&) const;

	extern template void System::SetPathVariable<complex_dbl>(complex_dbl const&) const;
	extern template void System::SetPathVariable<complex_mp>(complex_mp const&) const;

	extern template void System::SetAndReset<complex_dbl>(Vec<complex_dbl> const&, complex_dbl const&) const;
	extern template void System::SetAndReset<complex_mp>(Vec<complex_mp> const&, complex_mp const&) const;

	extern template void System::SetAndReset<complex_dbl>(Vec<complex_dbl> const&) const;
	extern template void System::SetAndReset<complex_mp>(Vec<complex_mp> const&) const;
	/// \endcond

}









#endif // for the ifndef include guards



