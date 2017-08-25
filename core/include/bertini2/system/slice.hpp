//This file is part of Bertini 2.
//
//slice.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//slice.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with slice.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file bertini2/system/slice.hpp 

\brief Provides the bertini::LinearSlice class.
*/


#ifndef BERTINI_SLICE_HPP
#define BERTINI_SLICE_HPP

#include "bertini2/function_tree.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"


namespace bertini {


	/**
	\brief Base class for other Slices.  

	\see LinearSlice
	*/
	class Slice
	{

	};



	/**
	\brief Slice an affine or projective space with a LinearSlice today!
	*/
	class LinearSlice : public Slice
	{


		mutable std::tuple<Mat<dbl>, Mat<mpfr> > coefficients_working_;
		Mat< mpfr > coefficients_highest_precision_; ///< the highest-precision coefficients for the patch

		mutable std::tuple<Vec<dbl>, Vec<mpfr> > constants_working_;
		Vec< mpfr > constants_highest_precision_; ///< the highest-precision coefficients for the patch

		VariableGroup sliced_vars_;
		unsigned num_dims_sliced_;
		mutable unsigned precision_; ///< the current working precision of the patch.

		bool is_homogeneous_;
	public:

		/**
		Produce a random real slice on a variable group, slicing a given number of dimensions.
		*/
		static LinearSlice RandomReal(VariableGroup const& v, unsigned dim, bool homogeneous = false, bool orthogonal = true)
		{
			return Make(v, dim, homogeneous, orthogonal, bertini::RandomReal);
		}

		/**
		\brief Generate a random complex slice.
		*/
		static LinearSlice RandomComplex(VariableGroup const& v, unsigned dim, bool homogeneous = false, bool orthogonal = true)
		{
			return Make(v, dim, homogeneous, orthogonal, bertini::RandomComplex);
		}


		/**
		\brief Evaluate the function values of the LinearSlice, in-place
		*/
		template<typename NumT>
		void Eval(Vec<NumT> & result, Vec<NumT> const& x) const
		{
			result = std::get<Mat<NumT> >(coefficients_working_) * x;

			if (!is_homogeneous_)
				result += std::get<Vec<NumT> >(constants_working_);
		}

		/**
		\brief Evaluate the function values of the LinearSlice
		*/
		template<typename NumT>
		Vec<NumT> Eval(Vec<NumT> const& x) const
		{
			if (!is_homogeneous_)
				return std::get<Mat<NumT> >(coefficients_working_) * x + std::get<Vec<NumT> >(constants_working_);
			else
				return std::get<Mat<NumT> >(coefficients_working_) * x;
		}


		/**
		\brief Evaluate the Jacobian of the LinearSlice, in-place
		*/
		template<typename NumT>
		void Jacobian(Mat<NumT> & result, Mat<NumT> const& x) const
		{
			result = std::get<Mat<NumT> >(coefficients_working_);
		}

		/**
		\brief Evaluate the Jacobian of the LinearSlice
		*/
		template<typename NumT>
		Mat<NumT> Jacobian(Mat<NumT> const& x) const
		{
			return std::get<Mat<NumT> >(coefficients_working_);
		}




		/**
		\brief Get the current precision of the slice.

		\return The current precision, in digits.
		*/
		unsigned Precision() const
		{
			return precision_;
		}

		/**
		\brief Set the precision of the slice.
	
		Copies the slice coefficients into correct precision for subsequent precision.

		\param new_precision The precision to change to.
		*/
		void Precision(unsigned new_precision) const
		{
			if (new_precision > DoublePrecision())
			{
				Mat<mpfr>& coefficients_mpfr = std::get<Mat<mpfr> >(coefficients_working_);
				Vec<mpfr>& constants_mpfr = std::get<Vec<mpfr> >(constants_working_);
				for (unsigned ii = 0; ii < Dimension(); ++ii)
				{
					for (unsigned jj=0; jj<NumVariables(); ++jj)
					{
						coefficients_mpfr(ii,jj).precision(new_precision);
						if (new_precision>precision_)
							coefficients_mpfr(ii,jj) = coefficients_highest_precision_(ii,jj);
					}
					
					if (!is_homogeneous_)
					{
						constants_mpfr(ii).precision(new_precision);
						if (new_precision>precision_)
							constants_mpfr(ii) = constants_highest_precision_(ii);
					}
				}
			}
			precision_ = new_precision;
		}

		/**
		\brief Get the dimension of the slice
		*/
		unsigned Dimension() const
		{
			return num_dims_sliced_;
		}

		/**
		\brief Get the number of variables sliced.  
		*/
		unsigned NumVariables() const
		{
			return sliced_vars_.size();
		}



		/** 
		\brief the default constructor for linear slices. 

		Make an empty linear slice. 
		*/
		LinearSlice() : 
			sliced_vars_(), 
			precision_(DefaultPrecision()), 
			num_dims_sliced_(0), 
			coefficients_highest_precision_(0, 0), 
			is_homogeneous_(false), 
			constants_highest_precision_(static_cast<unsigned>(0))
		{ 
			std::get<Mat<dbl> > (coefficients_working_).resize(Dimension(), NumVariables());
			std::get<Mat<mpfr> >(coefficients_working_).resize(Dimension(), NumVariables());

			std::get<Vec<dbl> > (constants_working_).resize(Dimension());
			std::get<Vec<mpfr> >(constants_working_).resize(Dimension());
		}

		/**
		\brief the constructor for linear slices.
		*/
		LinearSlice(VariableGroup const& v, unsigned dim, bool homogeneous) : sliced_vars_(v), precision_(DefaultPrecision()), num_dims_sliced_(dim), coefficients_highest_precision_(dim, v.size()), is_homogeneous_(homogeneous), constants_highest_precision_(dim)
		{ 
			std::get<Mat<dbl> > (coefficients_working_).resize(Dimension(), NumVariables());
			std::get<Mat<mpfr> >(coefficients_working_).resize(Dimension(), NumVariables());

			if (!homogeneous)
			{
				std::get<Vec<dbl> > (constants_working_).resize(Dimension());
				std::get<Vec<mpfr> >(constants_working_).resize(Dimension());
			}
		}


		/**
		\brief factory function for generating slices
		*/
		static
		LinearSlice Make(VariableGroup const& v, unsigned dim, bool homogeneous, bool orthogonal, std::function<void(mpfr&, unsigned)> gen)
		{
			LinearSlice s(v, dim, homogeneous);

			

			if (orthogonal)
			{	
				using std::min;
				using std::max;

				auto mindim = min(s.Dimension(),s.NumVariables());
				auto maxdim = max(s.Dimension(),s.NumVariables());

				bool need_transpose = s.Dimension() < s.NumVariables();

				s.coefficients_highest_precision_.resize(maxdim,mindim);

				for (unsigned ii(0); ii<maxdim; ++ii)
					for (unsigned jj(0); jj<mindim; ++jj)
						gen(s.coefficients_highest_precision_(ii,jj), MaxPrecisionAllowed());

				auto prev_precision = DefaultPrecision();
				DefaultPrecision(MaxPrecisionAllowed());

				auto QR_factorization = Eigen::HouseholderQR<Mat<mpfr> >(s.coefficients_highest_precision_);
				s.coefficients_highest_precision_ = QR_factorization.householderQ()*Mat<mpfr>::Identity(maxdim, mindim);
				
				if (need_transpose)
					s.coefficients_highest_precision_.transposeInPlace();

				DefaultPrecision(prev_precision);
			}
			else
			{
				for (unsigned ii(0); ii<s.Dimension(); ++ii)
					for (unsigned jj(0); jj<s.NumVariables(); ++jj)
						gen(s.coefficients_highest_precision_(ii,jj), MaxPrecisionAllowed());
			}

			for (unsigned ii(0); ii<s.Dimension(); ++ii)
				for (unsigned jj(0); jj<s.NumVariables(); ++jj)
				{
					std::get<Mat<dbl> >(s.coefficients_working_)(ii,jj) = dbl(s.coefficients_highest_precision_(ii,jj));
					std::get<Mat<mpfr> >(s.coefficients_working_)(ii,jj) = s.coefficients_highest_precision_(ii,jj);
				}

			if (!homogeneous)
			{
				s.constants_highest_precision_.resize(s.Dimension());
				std::get<Vec<dbl> >(s.constants_working_).resize(s.Dimension());
				std::get<Vec<mpfr> >(s.constants_working_).resize(s.Dimension());

				for (unsigned ii(0); ii<s.Dimension(); ++ii)
				{
					gen(s.constants_highest_precision_(ii), MaxPrecisionAllowed());
					std::get<Vec<dbl> >(s.constants_working_)(ii) = dbl(s.constants_highest_precision_(ii));
					std::get<Vec<mpfr> >(s.constants_working_)(ii) = s.constants_highest_precision_(ii);
				}
			}

			assert(s.coefficients_highest_precision_.rows()==s.Dimension());
			assert(s.coefficients_highest_precision_.cols()==s.NumVariables());

			if (!homogeneous)
				assert(s.constants_highest_precision_.size()==s.Dimension());
			else
				assert(s.constants_highest_precision_.size()==0);

			return s;
		}

	private:
		

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & precision_;

			ar & coefficients_highest_precision_;

			ar & std::get<0>(coefficients_working_);
			ar & std::get<1>(coefficients_working_);
			ar & sliced_vars_;
		}

		friend std::ostream& operator<<(std::ostream&, LinearSlice const&);
	};

	/**
	\brief Provides output streaming for LinearSlice
	*/
	std::ostream& operator<<(std::ostream& out, LinearSlice const& s);
} // re: namespace bertini 

#endif

