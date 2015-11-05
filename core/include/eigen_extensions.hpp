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

#include <eigen3/Eigen/Dense>


namespace bertini {

	using std::abs;
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

		return abs(numerator/denomenator) >= 1.0/Eigen::NumTraits<T>::dummy_precision();
	} 

	enum class MatrixSuccessCode
	{
		Success,
		LargeChange,
		SmallValue
	};

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


	template <typename NumberType>
	Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> KahanMatrix(unsigned int mat_size, NumberType c)
	{
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

}

namespace bertini{
	template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
	template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;
}

#endif

