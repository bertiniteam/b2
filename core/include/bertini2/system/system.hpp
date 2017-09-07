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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file system.hpp 

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
#include <boost/type_index.hpp>

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"


#include "bertini2/function_tree.hpp"
#include "bertini2/system/patch.hpp"

#include "bertini2/limbo.hpp"


#include <boost/archive/binary_oarchive.hpp>
	#include <boost/archive/binary_iarchive.hpp>
	#include <boost/iostreams/stream_buffer.hpp>
	#include <boost/iostreams/stream.hpp>
	#include <boost/iostreams/device/back_inserter.hpp>

namespace bertini {

	
	enum class JacobianEvalMethod
	{
		JacobianNode,
		Derivatives
		//StraightLineProgram not yet....
	};

	/**
	\brief Gets the default evaluation method for Jacobians.  One might be faster...
	*/
	JacobianEvalMethod DefaultJacobianEvalMethod();

	/**
	\brief Get the default value for whether a system should autosimplify.
	*/
	bool DefaultAutoSimplify();

	/**
	\brief The fundamental polynomial system class for Bertini2.
	
	 The fundamental polynomial system class for Bertini2.

	 Other System types are derived from this, but this class is not abstract.
	 */
	class System{
	
	public:
		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<node::Function>;
		using Var = std::shared_ptr<node::Variable>;
		using Nd = std::shared_ptr<node::Node>;
		using Jac = std::shared_ptr<node::Jacobian>;
		
		/**
		\brief The default constructor for a system.
		*/
		System() : is_differentiated_(false), have_path_variable_(false), have_ordering_(false), precision_(DefaultPrecision()), is_patched_(false)
		{}

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
		 \brief Evaluate the system using the previously set variable (and time) values, in place.

		It is up to YOU to ensure that the system's variables (and path variable) has been set prior to this function call.

		\return The function values of the system
		*/ 
		template<typename Derived>
		void EvalInPlace(Eigen::MatrixBase<Derived> & function_values) const
		{
			typedef typename Derived::Scalar T;

			if(function_values.size() < NumFunctions())
			{
				std::stringstream ss;
				ss << "trying to evaluate system in place, but number of input functions (" << function_values.size() << ") doesn't match number of system functions (" << NumFunctions() << ").";
				throw std::runtime_error(ss.str());
			}

			// the Reset() function call traverses the entire tree, resetting everything.
			// TODO: it has the unfortunate side effect of resetting constant functions, too.
			for (const auto& iter : functions_) 
				iter->Reset();


			unsigned counter(0);
			for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, counter++) {
				(*iter)->EvalInPlace<T>(function_values(counter));
			}

			if (IsPatched())
				patch_.EvalInPlace(function_values,
									std::get<Vec<T> >(current_variable_values_));
									// .segment(NumFunctions(),NumTotalVariableGroups())
			
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
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 \param variable_values The values of the variables, for the evaluation.
		 */
		template<typename Derived, typename OtherDerived>
		void EvalInPlace(Eigen::MatrixBase<Derived>& function_values, const Eigen::MatrixBase<OtherDerived>& variable_values) const
		{
			static_assert(std::is_same<typename Derived::Scalar,typename OtherDerived::Scalar>::value,"scalar types must match");

			if (variable_values.size()!=NumVariables())
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
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\param variable_values The values of the variables, for the evaluation.
		*/
		template<typename Derived>
		typename Derived::PlainObject Eval(const Eigen::MatrixBase<Derived>& variable_values) const
		{
			typedef typename Derived::Scalar T;

			if (variable_values.size()!=NumVariables())
			{
				std::stringstream ss;
				ss << "trying to evaluate system, but number of input variables (" << variable_values.size() << ") doesn't match number of system variables (" << NumVariables() << ").";
				throw std::runtime_error(ss.str());
			}
			if (have_path_variable_)
				throw std::runtime_error("not using a time value for evaluation of system, but path variable IS defined.");

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
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \todo The Eval() function for systems has the unfortunate side effect of resetting constant functions.  Modify the System class so that only certain parts of the tree get reset.
		 */
		template<typename Derived, typename OtherDerived, typename T>
		void EvalInPlace(Eigen::MatrixBase<Derived> & function_values, const Eigen::MatrixBase<OtherDerived>& variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");
			static_assert(std::is_same<typename OtherDerived::Scalar, T>::value, "scalar types must be the same");

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");

			SetVariables(variable_values.eval());//TODO: remove this eval
			SetPathVariable(path_variable_value);

			EvalInPlace(function_values);
		}

		
		
		
		




		/**
		 Evaluate the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \todo The Eval() function for systems has the unfortunate side effect of resetting constant functions.  Modify the System class so that only certain parts of the tree get reset.
		 */
		template<typename Derived, typename T>
		Vec<T> Eval(const Eigen::MatrixBase<Derived>& variable_values, const T & path_variable_value) const
		{

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
			if (!have_path_variable_)
				throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");

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

		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		*/
		template<typename Derived>
		void JacobianInPlace(Eigen::MatrixBase<Derived> & J) const
		{
			typedef typename Derived::Scalar T;

			if(J.rows() != NumTotalFunctions() || J.cols() != NumVariables())
			{
				throw std::runtime_error("trying to evaluate jacobian of system in place, but input J doesn't have right number of columns or rows");
			}
			
			const auto& vars = Variables();

			if (!is_differentiated_)
				Differentiate();
			
			

			switch (jacobian_eval_method_)
			{
				case JacobianEvalMethod::JacobianNode:
				{
					for (const auto& iter : jacobian_) 
						iter->Reset();

					for (int ii = 0; ii < NumFunctions(); ++ii)
						for (int jj = 0; jj < NumVariables(); ++jj)
							jacobian_[ii]->EvalJInPlace<T>(J(ii,jj),vars[jj]);
					break;
				}
				case JacobianEvalMethod::Derivatives:
				{
					for (const auto& iter : space_derivatives_) 
						iter->Reset();

					for (int jj = 0; jj < NumVariables(); ++jj)
						for (int ii = 0; ii < NumFunctions(); ++ii)
							space_derivatives_[ii+jj*NumFunctions()]->EvalInPlace<T>(J(ii,jj));
					break;
				}
			}
			
			if (IsPatched())
				patch_.JacobianInPlace(J,std::get<Vec<T> >(current_variable_values_));
			
		}

		
		
		
		
		


		/**
		Evaluate the Jacobian matrix of the system, using the previous space and time values.

		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
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
		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 */
		template<typename Derived, typename T>
		void JacobianInPlace(Eigen::MatrixBase<Derived> & J, const Vec<T> &  variable_values) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must match");

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");
			
			if (HavePathVariable())
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");
			
			SetVariables(variable_values);
			
			JacobianInPlace(J);
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

			if (HavePathVariable())
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J,variable_values);
			return J;
		}


		
		
		/**
		 Evaluate the Jacobian of the system, provided a path variable is defined for the system, in place.
		 
		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \return The Jacobian matrix.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		 */
		template<typename Derived, typename OtherDerived, typename T>
		void JacobianInPlace(Eigen::MatrixBase<Derived> & J, const Eigen::MatrixBase<OtherDerived> & variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");
			static_assert(std::is_same<typename OtherDerived::Scalar, T>::value, "scalar types must be the same");

			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");
			
			if (!HavePathVariable())
				throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");
			
			SetVariables(variable_values.eval()); // TODO: remove this eval
			SetPathVariable(path_variable_value);

			JacobianInPlace(J);
		}

		
		


		/**
		 Evaluate the Jacobian of the system, provided a path variable is defined for the system.

		 \throws std::runtime_error, if a path variable is NOT defined, and you passed it a value.  Also throws if the number of variables doesn't match.
		 \return The Jacobian matrix.
		 
		 \param variable_values The values of the variables, for the evaluation.
		 \param path_variable_value The current value of the path variable.

		 \tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
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
			if (variable_values.size()!=NumVariables())
				throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

			if (!HavePathVariable())
				throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");

			Mat<T> J(NumTotalFunctions(), NumVariables());
			JacobianInPlace(J,variable_values,path_variable_value);
			return J;
		}

		
		/**
		\brief Compute the time-derivative is a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename Derived, typename OtherDerived, typename T>
		void TimeDerivativeInPlace(Eigen::MatrixBase<Derived> & ds_dt, 
							const Eigen::MatrixBase<OtherDerived> & variable_values, 
							const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");
			static_assert(std::is_same<typename OtherDerived::Scalar, T>::value, "scalar types must be the same");

			if(ds_dt.size() < NumFunctions())
		{
				std::stringstream ss;
				ss << "trying to evaluate system in place, but number of input functions (" << ds_dt.size() << ") doesn't match number of system functions (" << NumFunctions() << ").";
				throw std::runtime_error(ss.str());
			}
			if (!HavePathVariable())
				throw std::runtime_error("computing time derivative of system with no path variable defined");

			if (!is_differentiated_)
				Differentiate();

			SetVariables(variable_values.eval()); //TODO: remove this eval()
			SetPathVariable(path_variable_value);

			switch (jacobian_eval_method_)
			{
				case JacobianEvalMethod::JacobianNode:
				{
					for (int ii = 0; ii < NumFunctions(); ++ii)
						jacobian_[ii]->Reset();

					for (int ii = 0; ii < NumFunctions(); ++ii)
						jacobian_[ii]->EvalJInPlace<T>(ds_dt(ii), path_variable_);
					break;
				}
				case JacobianEvalMethod::Derivatives:
				{
					for (int ii = 0; ii < NumFunctions(); ++ii)
						time_derivatives_[ii]->Reset();

					for (int ii = 0; ii < NumFunctions(); ++ii)
						time_derivatives_[ii]->EvalInPlace<T>(ds_dt(ii));
					break;
				}
			}

			// the patch doesn't move with time.  derivatives 0.
			if (IsPatched())
				for (int ii = 0; ii < NumTotalVariableGroups(); ++ii)
					ds_dt(ii+NumFunctions()) = T(0);
			
		}

		
		
		
		
		
		/**
		\brief Compute the time-derivative is a system. 
		
		If \f$S\f$ is the system, and \f$t\f$ is the path variable this computes \f$\frac{dS}{dt}\f$.

		\tparam T The number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.
		\throws std::runtime error if the system does not have a path variable defined.
		*/
		template<typename Derived, typename T>
		Vec<T> TimeDerivative(const Eigen::MatrixBase<Derived> & variable_values, const T & path_variable_value) const
		{
			static_assert(std::is_same<typename Derived::Scalar, T>::value, "scalar types must be the same");

			if (!HavePathVariable())
				throw std::runtime_error("computing time derivative of system with no path variable defined");


			Vec<T> ds_dt(NumTotalFunctions());
			TimeDerivativeInPlace(ds_dt, variable_values, path_variable_value);
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

			const auto& vars = Variables();

			#ifndef BERTINI_DISABLE_PRECISION_CHECKS
				if (!std::is_same<T,dbl>::value && (Precision(new_values) != this->precision()))
					throw std::runtime_error("precision of input point in SetVariables (" + std::to_string(Precision(new_values)) + ") must match the precision of the system (" + std::to_string(this->precision()) + ").");

				if (!std::is_same<T,dbl>::value && (vars[0]->node::NamedSymbol::precision() != this->precision()) )
					throw std::runtime_error("internally, precision of variables (" + std::to_string(vars[0]->node::NamedSymbol::precision()) + ") in SetVariables must match the precision of the system (" + std::to_string(this->precision()) + ").");
			#endif

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
		void SetPathVariable(T const& new_value) const
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

				return DehomogenizePointFIFO(x);
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
		\brief One of a family of functions indicating whether we can assume the system will always have uniform precision.

		\see PleaseAssumeUniformPrecision AssumeUniformPrecision IsAssumingUniformPrecision
		*/
		void DontAssumeUniformPrecision()
		{
			AssumeUniformPrecision(false);
		}

		void PleaseAssumeUniformPrecision()
		{
			AssumeUniformPrecision(true);
		}

		void AssumeUniformPrecision(bool val)
		{
			assume_uniform_precision_ = false;
		}

		/** 
		\brief yon getter for the obvious thing it gets
		*/
		auto IsAssumingUniformPrecision() const
		{
			return assume_uniform_precision_;
		}

		inline
		void PleaseAutoSimplify()
		{
			SetAutoSimplify(true);
		}

		inline
		void DontAutoSimplify()
		{
			SetAutoSimplify(false);
		}

		void SetAutoSimplify(bool val)
		{
			auto_simplify_ = val;
		}


		/**
		\brief Query the state of autosimplification
		*/
		auto IsAutoSimplifying() const
		{
			return auto_simplify_;
		}

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
	
		\tparam T the number-type for return.  Probably dbl=std::complex<double>, or mpfr=bertini::complex.

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
			unsigned ungrouped_variable_counter = 0;

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
		 Puts together the ordering of variables, and stores it internally.
		*/
		void ConstructOrdering() const;


		VariableGroup ungrouped_variables_; ///< ungrouped variable nodes.  Not in an affine variable group, not in a projective group.  Just hanging out, being a variable.
		std::vector< VariableGroup > variable_groups_; ///< Affine variable groups.  When system is homogenized, will have a corresponding homogenizing variable.
		std::vector< VariableGroup > hom_variable_groups_; ///< Homogeneous or projective variable groups.  System SHOULD be homogeneous with respect to these.  

		VariableGroup homogenizing_variables_; ///< homogenizing variables for the variable_groups.  


		bool have_path_variable_ = false; ///< Whether we have the variable or not.
		Var path_variable_; ///< the single path variable for this system.  Sometimes called time.
		
		VariableGroup implicit_parameters_; ///< Implicit parameters.  These don't depend on anything, and will be moved from one parameter point to another by the tracker.  They should be algebraically constrained by some equations.
		std::vector< Fn > explicit_parameters_; ///< Explicit parameters.  These should be functions of the path variable only, NOT of other variables.  

		std::vector< Fn > constant_subfunctions_; ///< degree-0 functions, depending on neither variables nor the path variable.
		std::vector< Fn > subfunctions_; ///< Any declared subfunctions for the system.  Can use these to ensure that complicated repeated structures are only created and evaluated once.
		std::vector< Fn > functions_; ///< The system's functions.
		
		class Patch patch_; ///< Patch on the variable groups.  Assumed to be in the same order as the time_order_of_variable_groups_ if the system uses FIFO ordering, or in same order as the AffHomUng variable groups if that is set.
		bool is_patched_ = false;	///< Indicator of whether the system has been patched.

		mutable std::vector< Jac > jacobian_; ///< The generated functions from differentiation.  Created when first call for a Jacobian matrix evaluation.

		mutable std::vector< Nd > space_derivatives_; ///< The generated functions from differentiation with respect to space.  in column-major order to be consistent with Eigen default order.  Created when first call for a Jacobian matrix evaluation.

		mutable std::vector< Nd > time_derivatives_; ///< The generated functions from differentiation with respect to time.  in column-major order to be consistent with Eigen default order.  Created when first call for a Jacobian matrix evaluation.

		mutable bool is_differentiated_ = false; ///< indicator for whether the jacobian tree has been populated.


		std::vector< VariableGroupType > time_order_of_variable_groups_;

		mutable std::tuple< Vec<dbl>, Vec<mpfr> > current_variable_values_;

		mutable VariableGroup variable_ordering_; ///< The assembled ordering of the variables in the system.
		mutable bool have_ordering_ = false;

		mutable unsigned precision_; ///< the current working precision of the system 

		bool assume_uniform_precision_ = false; ///< a bit, setting whether we can assume the system is in uniform precision.  if you are doing things that will allow pieces of the system to drift in terms of precision, then you should not assume this.  \see AssumeUniformPrecision

		JacobianEvalMethod jacobian_eval_method_ = DefaultJacobianEvalMethod(); ///< an enum class value, indicating which method of evaluation should be used.

		bool auto_simplify_ = DefaultAutoSimplify();

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

			ar & variable_ordering_;
			ar & have_ordering_;

			ar & implicit_parameters_;
			ar & explicit_parameters_;

			ar & constant_subfunctions_;
			ar & subfunctions_;
			ar & functions_;

			ar & is_differentiated_;
			ar & jacobian_;

			ar & space_derivatives_;
			ar & time_derivatives_;
			
			ar & precision_;
			ar & is_patched_;
			ar & patch_;

			ar & assume_uniform_precision_;
			ar & jacobian_eval_method_;
		}

	};


	/**
	\brief Contcatenate two compatible systems.

	Two systems are compatible for concatenation if they have the same variable structure, and if they have the same patch (if patched).

	\param sys1 The top system.
	\param sys2 The bottom system.

	If both patched both must have same patch.  If not both are patched, then the patch will propagate to the returned system. 

	If the two patches have differing variable orderings, the call to Concatenate will throw.
	*/
	System Concatenate(System sys1, System const& sys2);
	


	/**
	\brief Do a deep clone of the system.  This includes the entire structure, variables, etc.  everything.
	*/
	System Clone(System const& sys);

	/**
	\brief Free form function for simplifying systems.
	*/
	void Simplify(System & sys);
}









#endif // for the ifndef include guards



