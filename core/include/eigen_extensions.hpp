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

#ifndef BERTINI_EIGEN_EXTENSIONS
#define BERTINI_EIGEN_EXTENSIONS

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


}



#endif

