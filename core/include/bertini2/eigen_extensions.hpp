//This file is part of Bertini 2.0.
//
//eigen_extensions.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//eigen_extensions.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with .  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015
//
//  This file contains extensions to types provided by Eigen.
//
//

/**
\file eigen_extensions.hpp 

\brief Bertini extensions of the Eigen linear algebra library.
*/

#ifndef BERTINI_EIGEN_EXTENSIONS_HPP
#define BERTINI_EIGEN_EXTENSIONS_HPP


#include "num_traits.hpp"


#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>
// reopen the Eigen namespace to inject this struct.
namespace Eigen {
	
	using mpfr_float = bertini::mpfr_float;
	template<> struct NumTraits<mpfr_float> : GenericNumTraits<mpfr_float> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
		
		typedef mpfr_float Real;
		typedef mpfr_float NonInteger;
		typedef mpfr_float Nested;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 20,
			AddCost = 30,
			MulCost = 40
		};
		
		
		inline static Real highest() {
			
			return (mpfr_float(1) - epsilon()) * pow(mpfr_float(2),mpfr_get_emax()-1);//);//mpfr_float::default_precision());
		}
		
		inline static Real lowest() {
			return -highest();
		}
		
		inline static Real dummy_precision()
		{
			return pow( mpfr_float(10),-int(mpfr_float::default_precision()-3));
		}
		
		inline static Real epsilon()
		{
			return pow(mpfr_float(10),-int(mpfr_float::default_precision()));
		}
		//http://www.manpagez.com/info/mpfr/mpfr-2.3.2/mpfr_31.php
	};

}





// reopen the Eigen namespace to inject this struct.
namespace Eigen {
	
	
	
	/**
	 \brief This templated struct permits us to use the bertini::complex type in Eigen matrices.

	 Provides methods to get the epsilon, dummy_precision, lowest, highest functions, largely by inheritance from the NumTraits<mpfr_float> contained in mpfr_extensions.
	 */
	template<> struct NumTraits<bertini::complex> : NumTraits<bertini::mpfr_float> 
	{
		using mpfr_float = bertini::mpfr_float;

		typedef mpfr_float Real;
		typedef mpfr_float NonInteger;
		typedef bertini::complex Nested;// Nested;
		enum {
			IsComplex = 1,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 2 * NumTraits<Real>::ReadCost,
			AddCost = 2 * NumTraits<Real>::AddCost,
			MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
		};
		
	};

	namespace internal {
		template<>
		struct abs2_impl<bertini::complex>
		{
			static inline mpfr_float run(const bertini::complex& x)
			{
				return real(x)*real(x) + imag(x)*imag(x);
			}
		};


		template<> inline bertini::complex random<bertini::complex>()
		{
			return bertini::complex::rand();
		}

		template<> inline bertini::complex random<bertini::complex>(const bertini::complex& a, const bertini::complex& b)
		{
			return a + (b-a) * random<bertini::complex>();
		}

		template<>
		struct conj_helper<bertini::complex, bertini::complex, false, true>
		{
			typedef bertini::complex Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }

			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(numext::real(x)*numext::real(y) + numext::imag(x)*numext::imag(y), numext::imag(x)*numext::real(y) - numext::real(x)*numext::imag(y)); }
		};

		template<>
		struct conj_helper<bertini::complex, bertini::complex, true, false>
		{
			typedef bertini::complex Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }

			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(numext::real(x)*numext::real(y) + numext::imag(x)*numext::imag(y), numext::real(x)*numext::imag(y) - numext::imag(x)*numext::real(y)); }
		};

	} // re: namespace internal
} // re: namespace Eigen

#include <Eigen/Dense>

namespace bertini {

	template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
	template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;


	
	/**
	 Test a numbers being very small.

	 Compares the number against machine epsilon (or software epsilon if a multiple precision type), times 100.  

	 \note Machine epsilon for doubles is about 1e-16, for mpfr_float, it's 10^-current precision.
	 
	 \param testme The number to test

	 \return true, if the number is very small.  False otherwise.
	*/
	template<typename T>
	bool IsSmallValue(T const& testme)
	{
		using std::abs;
		return abs(testme) <= Eigen::NumTraits<T>::epsilon()*100;
	} 

	/**
	 Test two numbers for having large ratio.
	
	 The basis for comparison is Eigen's fuzzy precision, Eigen::NumTraits<T>::dummy_precision();
	 
	 \note For doubles, the threshold is 1e-12, for mpfr_float is 1e3*current precision.

	 \param numerator The numerator for the ratio.
	 \param denomenator The denomenator for the ratio.

	 \return true, if the ratio is very large, false otherwise
	*/
	template<typename T>
	bool IsLargeChange(T const& numerator ,T const& denomenator)
	{
		using std::abs;
		return abs(numerator/denomenator) >= 1.0/Eigen::NumTraits<T>::dummy_precision();
	} 

	enum class MatrixSuccessCode
	{
		Success,
		LargeChange,
		SmallValue
	};

	/**
	\brief Check the diagonal elements of an LU decomposition for small values and large ratios.  

	\return Success if things are ok.  LargeChange or SmallValue if one is found.

	This function requires a square non-empty matrix.

	\tparam Derived Matrix type from Eigen.
	*/
	template <typename Derived>
	MatrixSuccessCode LUPartialPivotDecompositionSuccessful(Eigen::MatrixBase<Derived> const& LU)
	{
		#ifndef BERTINI_DISABLE_ASSERTS
			assert(LU.rows()==LU.cols() && "non-square matrix in LUPartialPivotDecompositionSuccessful");
			assert(LU.rows()>0 && "empty matrix in LUPartialPivotDecompositionSuccessful");
		#endif

			// this loop won't test entry (0,0).  it's tested separately after.
		for (unsigned int ii = LU.rows()-1; ii > 0; ii--)
		{
			if (IsSmallValue(LU(ii,ii))) 
			{
				return MatrixSuccessCode::SmallValue;
			}

			if (IsLargeChange(LU(ii-1,ii-1),LU(ii,ii)))
			{
				return MatrixSuccessCode::LargeChange;
			}
		}

		// this line is the reason for the above assert on non-empty matrix.
		if (IsSmallValue(LU(0,0))) 
		{
			return MatrixSuccessCode::SmallValue;
		}

		return MatrixSuccessCode::Success;
	}

	/**
	\brief Make a Kahan matrix with a given number type.
	*/
	template <typename NumberType>
	Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> KahanMatrix(unsigned int mat_size, NumberType c)
	{
		using std::sqrt;
		NumberType s, scale(1.0);
		s = sqrt( (NumberType(1.0)-c) * (NumberType(1.0)+c) );
		
		
		Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size);
		
		
		for (unsigned int ii=0; ii<mat_size; ii++) {
			for (unsigned int jj=0; jj<ii; jj++) {
				A(ii,jj) = NumberType(0.0);
			}
			for (unsigned int jj=ii; jj<mat_size; jj++) {
				A(ii,jj) = -c * NumberType(1.0);
			}
		}
		
		
		for (unsigned int ii=0; ii<mat_size; ii++) {
			A(ii,ii) += NumberType(1)+c;
		}
		
		for (unsigned int jj=0; jj<mat_size; jj++){
			for (unsigned int kk=0;kk<mat_size;kk++){
				A(jj,kk) *= scale;
			}
			scale *= s;
		}
		
		for (unsigned int jj=0;jj<mat_size;jj++){
			for (unsigned int kk=0;kk<mat_size;kk++){
				A(kk,jj)/= NumberType(jj) + NumberType(1);
			}
		}
		return A;
	}


	/**
	\brief Make a random matrix of units (numbers with norm 1).

	\return The random matrix of units.

	\param rows The number of rows.
	\param cols The number of columns.
	\tparam NumberType the type of number to fill the matrix with.
	*/
	template <typename NumberType>
	Mat<NumberType> RandomOfUnits(uint rows, uint cols)
	{
		return Mat<NumberType>(rows,cols).unaryExpr([](NumberType const& x) { return RandomUnit<NumberType>(); });
	}

	/**
	\brief Make a random vector of units (numbers with norm 1).

	\return The random vector of units.

	\param size The length of the vector.
	\tparam NumberType the type of number to fill the vector with.
	*/
	template <typename NumberType>
	Vec<NumberType> RandomOfUnits(uint size)
	{
		return Vec<NumberType>(size).unaryExpr([](NumberType const& x) { return RandomUnit<NumberType>(); });
	}

}


// the following code comes from
// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
// and adds support for serialization of Eigen types
//
// question asked by user Gabriel and answered by Shmuel Levine
// answer code repeated here verbatim.
// please update this comment if this code is changed, 
// and post the modifications to the above referenced post on SO.
namespace boost{
    namespace serialization{

        template<   class Archive, 
                    class S, 
                    int Rows_, 
                    int Cols_, 
                    int Ops_, 
                    int MaxRows_, 
                    int MaxCols_>
        inline void save(
            Archive & ar, 
            const Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g, 
            const unsigned int version)
            {
                int rows = g.rows();
                int cols = g.cols();

                ar & rows;
                ar & cols;
                ar & boost::serialization::make_array(g.data(), rows * cols);
            }

        template<   class Archive, 
                    class S, 
                    int Rows_,
                    int Cols_,
                    int Ops_, 
                    int MaxRows_, 
                    int MaxCols_>
        inline void load(
            Archive & ar, 
            Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g, 
            const unsigned int version)
        {
            int rows, cols;
            ar & rows;
            ar & cols;
            g.resize(rows, cols);
            ar & boost::serialization::make_array(g.data(), rows * cols);
        }

        template<   class Archive, 
                    class S, 
                    int Rows_, 
                    int Cols_, 
                    int Ops_, 
                    int MaxRows_, 
                    int MaxCols_>
        inline void serialize(
            Archive & ar, 
            Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g, 
            const unsigned int version)
        {
            split_free(ar, g, version);
        }


    } // namespace serialization
} // namespace boost




#endif

