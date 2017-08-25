//This file is part of Bertini 2.
//
//include/bertini2/system/patch.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/system/patch.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/system/patch.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file include/bertini2/system/patch.hpp 

\brief Provides the bertini::Patch class.
*/

#ifndef BERTINI_PATCH_HPP
#define BERTINI_PATCH_HPP

#include "bertini2/mpfr_complex.hpp"

#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"

#include <vector>

#include <Eigen/Dense>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/vector.hpp>

namespace bertini {




	/**
	\class Patch

	\brief An affine chart on a product of projective spaces.

	This class provides an adjustable-precision patch on a working space for solving a polynomial system using homotopy continuation.

	The constant values for the patch are unity.

	## Example code, not using a bertini::System to generate the patch. 
	\code

	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<mpfr> v(5);
	v << mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1);

	auto v_rescaled = p.RescalePoint(v);

	auto f = p.Eval(v_rescaled);
	
	\endcode
	
	At the end of the above example, the patch evaluates to 0 in the container \f$f\f$.

	## Other ways of patching

	Patches can also be applied to a system, which should first be homogenized.
	*/
	class Patch
	{

	public:

		/**
		\brief Default constructor.

		Makes an empty patch.
		*/
		Patch() : precision_(DefaultPrecision())
		{}


		/**
		Custom copy constructor, ensuring the max-precision coefficients are copied in highest precision
		*/
		Patch(Patch const& other)
		{	
			variable_group_sizes_ = other.variable_group_sizes_;
			precision_ = DefaultPrecision();

			// a little shorthand unpacking the tuple
			std::vector<Vec<mpfr> >& coefficients_mpfr = std::get<std::vector<Vec<mpfr> > >(this->coefficients_working_);
			std::vector<Vec<dbl> >& coefficients_dbl = std::get<std::vector<Vec<dbl> > >(this->coefficients_working_);

			coefficients_highest_precision_.resize(other.NumVariableGroups());
			coefficients_mpfr.resize(variable_group_sizes_.size());
			coefficients_dbl.resize(variable_group_sizes_.size());

			for (unsigned ii(0); ii<other.NumVariableGroups(); ++ii)
			{
				auto curr_size = variable_group_sizes_[ii];

				coefficients_highest_precision_[ii].resize(curr_size);
				coefficients_dbl[ii].resize(curr_size);
				coefficients_mpfr[ii].resize(curr_size);

				for (unsigned jj(0); jj<curr_size; jj++)
				{
					coefficients_highest_precision_[ii](jj).precision(other.coefficients_highest_precision_[ii](jj).precision());

					coefficients_highest_precision_[ii](jj) = other.coefficients_highest_precision_[ii](jj);

					coefficients_dbl[ii](jj) = dbl(coefficients_highest_precision_[ii](jj));
					coefficients_mpfr[ii](jj) = mpfr(coefficients_highest_precision_[ii](jj));

					assert(coefficients_highest_precision_[ii](jj) == other.coefficients_highest_precision_[ii](jj));
				}
			}
		}


		/**
		\brief Constructor making a random complex patch on a space whose structure is described by the input argument.

		The sizes input give the number of total variables, including homogenizing variables, for the product of spaces forming the total space to be patched.  The patch has no idea whether the underlying space is projective or affine -- that is handled somewhere else, likely in a bertini::System.
		
		The initial precision of the patch is set to current default precision, and the precision of the highest precision coefficients are set to the current default as well.

		\param sizes The sizes of the variable groups, including homogenizing variables if present.
		*/
		Patch(std::vector<unsigned> const& sizes) : variable_group_sizes_(sizes), coefficients_highest_precision_(sizes.size()), precision_(DefaultPrecision())
		{
			using bertini::Precision;
			using bertini::RandomComplex;

			std::vector<Vec<mpfr> >& coefficients_mpfr = std::get<std::vector<Vec<mpfr> > >(coefficients_working_);
			std::vector<Vec<dbl> >& coefficients_dbl = std::get<std::vector<Vec<dbl> > >(coefficients_working_);

			coefficients_highest_precision_.resize(sizes.size());
			coefficients_dbl.resize(sizes.size());
			coefficients_mpfr.resize(sizes.size());

			for (size_t ii=0; ii<sizes.size(); ++ii)
			{
				coefficients_highest_precision_[ii].resize(sizes[ii]);
				coefficients_dbl[ii].resize(sizes[ii]);

				for (unsigned jj=0; jj<sizes[ii]; ++jj)
				{
					RandomComplex(coefficients_highest_precision_[ii](jj), MaxPrecisionAllowed());
					coefficients_dbl[ii](jj) = dbl(coefficients_highest_precision_[ii](jj));
				}

				coefficients_mpfr[ii] = coefficients_highest_precision_[ii]; // copy into current default precision
				Precision(coefficients_mpfr[ii],precision_);

				assert(Precision(coefficients_mpfr[ii](0))==precision_);
			}
		}


		/**
		\brief Construct a random complex patch on a space
		*/
		static Patch Random(std::vector<unsigned> const& sizes)
		{
			return Patch(sizes);
		}


		/**
		\brief Construct random REAL patch on a space.
		*/
		static Patch RandomReal(std::vector<unsigned> const& sizes)
		{	
			using bertini::RandomReal;
			using bertini::Precision;

			Patch p;

			p.variable_group_sizes_ = sizes;

			

			std::vector<Vec<mpfr> >& coefficients_mpfr = std::get<std::vector<Vec<mpfr> > >(p.coefficients_working_);
			std::vector<Vec<dbl> >& coefficients_dbl = std::get<std::vector<Vec<dbl> > >(p.coefficients_working_);

			p.coefficients_highest_precision_.resize(sizes.size());
			coefficients_mpfr.resize(sizes.size());
			coefficients_dbl.resize(sizes.size());
			
			for (size_t ii=0; ii<sizes.size(); ++ii)
			{
				p.coefficients_highest_precision_[ii].resize(sizes[ii]);
				for (unsigned jj=0; jj<sizes[ii]; ++jj)
					RandomReal(p.coefficients_highest_precision_[ii](jj), MaxPrecisionAllowed());

				coefficients_mpfr[ii] = p.coefficients_highest_precision_[ii]; 
				Precision(coefficients_mpfr[ii],DefaultPrecision());
				assert(Precision(coefficients_mpfr[ii](0))==DefaultPrecision());

				coefficients_dbl[ii].resize(sizes[ii]);
				for (unsigned jj=0; jj<sizes[ii]; ++jj)
					coefficients_dbl[ii](jj) = dbl(p.coefficients_highest_precision_[ii](jj));
			}

			return p;
		}


		/**
		\brief Get the current precision of the patch.

		\return The current precision, in digits.
		*/
		unsigned Precision() const
		{
			return precision_;
		}

		/**
		\brief Set the precision of the patch.
	
		Copies the patch coefficients into correct precision for subsequent precision.

		\param new_precision The precision to change to.
		*/
		void Precision(unsigned new_precision) const
		{
			using bertini::Precision;
			std::vector<Vec<mpfr> >& coefficients_mpfr = std::get<std::vector<Vec<mpfr> > >(coefficients_working_);

			for (unsigned ii = 0; ii < NumVariableGroups(); ++ii)
			{
				for (unsigned jj=0; jj<variable_group_sizes_[ii]; ++jj)
				{
					if (new_precision>precision_)
					{
						coefficients_mpfr[ii](jj) = coefficients_highest_precision_[ii](jj);
						coefficients_mpfr[ii](jj).precision(new_precision);
					}
					else
						coefficients_mpfr[ii](jj).precision(new_precision);

					assert(Precision(coefficients_mpfr[ii](jj))==new_precision);
				}
			}
			
			precision_ = new_precision;
		}


		

		

		/**
		\brief Evaluate the patch at a point, in place.

		\param function_values The vector to populate.  Must be at least as long as the number of variable groups.
		\param x The current space point at which to evaluate.

		\todo Rewrite this code to use Eigen sub-vectors, if possible.  If not, take this off the todo list.  See 
		http://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html

		*/
		template<typename Derived, typename T>
		void EvalInPlace(Eigen::MatrixBase<Derived> & function_values, Vec<T> const& x) const
		{
			static_assert(std::is_same<typename Derived::Scalar,T>::value,"scalar types must match");

			#ifndef BERTINI_DISABLE_ASSERTS
			assert(function_values.size()>=NumVariableGroups() && "function values must be of length at least as long as the number of variable groups");
//			assert((bertini::Precision(x(0))==DoublePrecision() || bertini::Precision(x(0)) == Precision())
//			 		&& "precision of input vector must match current working precision of patch during evaluation"
//			 	  );
			#endif

			// unpack from the tuple of working coefficients
			const std::vector<Vec<T> >& coefficients = std::get<std::vector<Vec<T> > >(coefficients_working_);

			unsigned offset(function_values.size() - NumVariableGroups()); // by precondition this number is at least 0.  the precondition is ensured by the public wrapper
			unsigned counter(0);
			for (unsigned ii = 0; ii < NumVariableGroups(); ++ii)
			{
				T& value = function_values(ii+offset);
				value = T(-1);
				for (unsigned jj=0; jj<variable_group_sizes_[ii]; ++jj)
				{	
					value += x(counter)*coefficients[ii](jj);
					counter++;
				}
			}
		}

		/**
		\brief Evaluate the patch at a point.

		\return The values of the patch at the point.
		\param x The current space point at which to evaluate.
		*/
		template<typename T>
		Vec<T> Eval(Vec<T> const& x) const
		{
			Vec<T> function_values = Vec<T>::Zero(NumVariableGroups());
			EvalInPlace(function_values, x);
			return function_values;
		}

		/**
		\brief Evaluate the Jacobian matrix, in place.

		\param jacobian Matrix to populate with the Jacobian.  Must be large enough (NumVariableGroups x NumVariables).
		\param x Point at which to evaluate.  Not technically needed, because the Jacobian is simply the matrix of coefficients.

		\todo Rewrite this code to use Eigen sub-vectors, if possible.  If not, take this off the todo list.  See 
		http://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html

		*/
		template<typename Derived, typename T>
		void JacobianInPlace(Eigen::MatrixBase<Derived> & jacobian, Vec<T> const& x) const
		{
			static_assert(std::is_same<typename Derived::Scalar,T>::value,"scalar types must match");


			#ifndef BERTINI_DISABLE_ASSERTS
			assert(jacobian.rows()>=NumVariableGroups() && "input jacobian must have at least as many rows as variable groups");
			assert(jacobian.cols()==NumVariables() && "input jacobian must have as many columns as the patch has variables");
			assert(
			       (bertini::Precision(x(0))==DoublePrecision() || bertini::Precision(x(0)) == Precision())  
			       	    && "precision of input vector must match current working precision of patch during evaluation"
			       );
			#endif
			
			const std::vector<Vec<T> >& coefficients = std::get<std::vector<Vec<T> > >(coefficients_working_);

			unsigned offset(jacobian.rows() - NumVariableGroups()); // by precondition this number is at least 0.  the precondition is ensured by the public wrapper
			unsigned counter(0);
			for (unsigned ii = 0; ii < NumVariableGroups(); ++ii)
				for (unsigned jj=0; jj<variable_group_sizes_[ii]; ++jj)
					jacobian(ii+offset,counter++) = coefficients[ii](jj);
		}

		/**
		\brief Evaluate the Jacobian matrix, in place.

		\param x Point at which to evaluate.  Not technically needed, because the Jacobian is simply the matrix of coefficients.
		\return Jacobian matrix, of size (NumVariableGroups x NumVariables).
		*/
		template<typename T>
		Mat<T> Jacobian(Vec<T> const& x) const
		{
			Mat<T> jacobian = Mat<T>::Zero(NumVariableGroups(), NumVariables());
			JacobianInPlace(jacobian, x);
			return jacobian;
		}

		
		/**
		\brief Rescale a point so that it satisfies the patch equations herein contained.

		The point must be as long as there are total variables.
		
		This function modifies the vector in place, so ensure you have a copy stored somewhere else if you feel it necessary to still have the pre-rescaled vector around.

		\param x The point to rescale
		\tparam T The number type.
		*/
		template<typename T>
		void RescalePointToFitInPlace(Vec<T> & x) const
		{
			#ifndef BERTINI_DISABLE_ASSERTS
				assert(x.size() == NumVariables() && "input point for rescaling to fit a patch must have same length as total number of variables being patched, in all variable groups.");
				assert((bertini::Precision(x(0))==DoublePrecision() || bertini::Precision(x(0)) == Precision())
						&& "precision of input vector must match current working precision of patch during rescaling"
					   );
			#endif

			const std::vector<Vec<T> >& coefficients = std::get<std::vector<Vec<T> > >(coefficients_working_);

			unsigned starting_index_counter(0);
			for (unsigned ii=0; ii<NumVariableGroups(); ii++)
			{
				auto subvec = x.segment(starting_index_counter,variable_group_sizes_[ii]);
				subvec /= subvec.transpose() * coefficients[ii];
				starting_index_counter += variable_group_sizes_[ii];
			}
		}

		/**
		\brief Rescale a point so that it satisfies the patch equations herein contained.

		\return Point \f$x\f$ rescaled to fit the patch.  \f$x\f$ cannot be the zero vector -- zero is a degenerate point in projective space.
		\param x The point to rescale.
		\tparam T The number type.
		*/
		template <typename T>
		Vec<T> RescalePoint(Vec<T> const& x) const
		{
			Vec<T> x_rescaled = x;
			RescalePointToFitInPlace(x_rescaled);
			return x_rescaled;
		}

		/**
		\brief Get the number of variable groups

		\return The number of variable groups in the problem.
		*/
		unsigned NumVariableGroups() const
		{
			return variable_group_sizes_.size();
		}


		/**
		\brief Get the number of variables in the patch.  
		*/
		unsigned NumVariables() const
		{
			unsigned num_vars(0);
			for (auto v : variable_group_sizes_)
				num_vars += v;
			return num_vars;
		}


		friend std::ostream& operator<<(std::ostream & out, Patch const& p)
		{
			out << p.NumVariableGroups() << " variable groups being patched\n";
			unsigned counter(0);
			for (auto& c : p.coefficients_highest_precision_)
			{
				out << "patch " << counter++ << " has " << c.size() << " coefficients:\n";
				out << c << "\n";
			}
			out << "current patch precision: " << p.Precision() << "\n";
			return out;
		}

		/**
		\brief Check whether two patches are the same. 

		\return true If have the same variable structure and coefficients.  Otherwise, false.
		*/
		bool operator==(Patch const rhs) const
		{
			if (NumVariableGroups()!=rhs.NumVariableGroups())
				return false;
			if (NumVariables()!=rhs.NumVariables())
				return false;
			for (unsigned ii=0; ii<NumVariableGroups(); ii++)
				if (variable_group_sizes_[ii]!=rhs.variable_group_sizes_[ii])
					return false;
			for (unsigned ii=0; ii<NumVariableGroups(); ii++)
				if (coefficients_highest_precision_[ii]!=rhs.coefficients_highest_precision_[ii])
					return false;

			return true;
		}

		/**
		\brief Check whether two patches are different.
		*/
		bool operator!=(Patch const& rhs) const
		{
			return !(*this==rhs);
		}

	private:

		/////////////////
		//
		//    Data members
		//
		//////////////////

		std::vector< Vec< mpfr > > coefficients_highest_precision_; ///< the highest-precision coefficients for the patch

		mutable std::tuple< std::vector< Vec< mpfr > >, std::vector< Vec< dbl > > > coefficients_working_; ///< the current working coefficients of the patch.  changing precision affects these, particularly the mpfr coefficients, which are down-sampled from the highest_precision coefficients.  the doubles are only down-sampled at time of creation or modification.

		std::vector<unsigned> variable_group_sizes_; ///< the sizes of the groups.  In principle, these must be at least 2.

		mutable unsigned precision_; ///< the current working precision of the patch.

		// add serialization support through boost.

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & precision_;

			ar & coefficients_highest_precision_;

			ar & std::get<0>(coefficients_working_);
			ar & std::get<1>(coefficients_working_);
			ar & variable_group_sizes_;
			
		}

	};


}


#endif  // re: include guards


