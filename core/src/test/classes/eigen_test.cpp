//This file is part of Bertini 2.0.
//
//eigen_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//eigen_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with eigen_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015




#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <boost/test/unit_test.hpp>


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>

#include "limbo.hpp"

#include "eigen_extensions.hpp"

#include "mpfr_complex.hpp"




extern double relaxed_threshold_clearance_d;
extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned CLASS_TEST_MPFR_DEFAULT_DIGITS;


using mpfr_float = boost::multiprecision::mpfr_float;

using bertini::KahanMatrix;







BOOST_AUTO_TEST_SUITE(kahan_matrix_solving_LU)



	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_double) {
		unsigned int size = 10;
		srand(2);  rand();
		
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A = KahanMatrix(size, 0.285), B(size,size), C;
		
		for (int ii=0; ii<size; ii++)
			for (int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -1.0/(ii+1) + double(rand()) /  RAND_MAX : B(ii,jj) = 0;
		
		
//		boost::timer::auto_cpu_timer t;
		C = A.lu().solve(B);
//		std::cout << C << std::endl;
//		std::cout << "pure double time to solve:" << std::endl;
		//add statement on the value of C to actually test

	}



	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_16)
	{
		
		unsigned int size = 10;
		boost::multiprecision::mpfr_float::default_precision(16);
		
		srand(2);  rand();
		
		Eigen::Matrix<boost::multiprecision::mpfr_float, Eigen::Dynamic, Eigen::Dynamic> A =
			KahanMatrix(size, boost::multiprecision::mpfr_float(0.285)), B(size,size), C;
		
		
		for (int ii=0; ii<size; ii++)
			for (int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -1.0/(ii+1.0) + boost::multiprecision::mpfr_float(rand()) /  boost::multiprecision::mpfr_float(RAND_MAX) : B(ii,jj) = 0;
		
		
//		boost::timer::auto_cpu_timer t;
		C = A.lu().solve(B);
//		std::cout << C << std::endl;
//		std::cout << "mpfr_float pure 16 time to solve:" << std::endl;
	}




	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_100)
	{
		unsigned int size = 10;
		srand(2);  rand();
		
		using mpfr = boost::multiprecision::mpfr_float;
		using mpfr_matrix = Eigen::Matrix<mpfr, Eigen::Dynamic, Eigen::Dynamic>;

		mpfr::default_precision(50);
		
		mpfr_matrix A = KahanMatrix(size, mpfr(0.285)), B(size,size), C;
		
		for (int ii=0; ii<size; ii++){
			for (int jj=0; jj<size; jj++){
				(jj!=ii) ? B(ii,jj) = -1.0/(ii+1) + mpfr(rand()) /  mpfr(RAND_MAX) : B(ii,jj) = mpfr(0.0);
			}
		}
		
		
	
		
//		boost::timer::auto_cpu_timer t;
		C = A.lu().solve(B);
//		std::cout << C << std::endl;
//		std::cout << "mpfr_float pure 100 time to solve:" << std::endl;
		
	}


	
	

	
	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_standardcomplex)
	{
		unsigned int size = 10;
		srand(2);  rand();
		
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A =
		KahanMatrix(size, std::complex<double>(0.285)), B(size,size), C;
		
		for (int ii=0; ii<size; ii++)
			for (int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -1.0/(ii+1) + double(rand()) / double(RAND_MAX) : B(ii,jj) = 0;
		
		
		
		
//		boost::timer::auto_cpu_timer t;
		C = A.lu().solve(B);
//		std::cout << C << std::endl;
//		std::cout << "std::complex time to solve:" << std::endl;
	}
	
	
	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_bertinicomplex_100)
	{
		
		unsigned int size = 10;
		boost::multiprecision::mpfr_float::default_precision(100);
		
		srand(2);  rand();
		
		Eigen::Matrix<bertini::complex, Eigen::Dynamic, Eigen::Dynamic> A =
		KahanMatrix(size, bertini::complex("0.285","0.0")), B(size,size), C;
		
		for (int ii=0; ii<size; ii++)
			for (int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = bertini::complex( bertini::complex(-1)/bertini::complex(ii+1) + bertini::complex(rand()) / bertini::complex(RAND_MAX)) : B(ii,jj) = bertini::complex(0);
		
//		boost::timer::auto_cpu_timer t;
		C = A.lu().solve(B);
		
		
//		std::cout << C << std::endl;
//		std::cout << "bertini::complex precision 100 time to solve:" << std::endl;
	}

		
		
		
		
//		BOOST_AUTO_TEST_CASE(mpfr_float_num_traits){
//			
//			std::cout << Eigen::NumTraits<boost::multiprecision::mpfr_float>::highest() << std::endl;
//			std::cout << Eigen::NumTraits<boost::multiprecision::mpfr_float>::lowest() << std::endl;
//			std::cout << Eigen::NumTraits<boost::multiprecision::mpfr_float>::dummy_precision() << std::endl;
//			std::cout << Eigen::NumTraits<boost::multiprecision::mpfr_float>::epsilon() << std::endl;
//			
//		}
	

	BOOST_AUTO_TEST_CASE(eigen_partial_pivot_solve_singular_matrix)
	{

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(2,2), B(2,1);

		A << 1, 1, 0, 0;

		B << 0.5, 1;

		auto LU = A.lu();


		auto C = LU.solve(B);

	}


	BOOST_AUTO_TEST_CASE(eigen_partial_pivot_solve_near_singular_matrix_double)
	{

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(2,2), B(2,1);

		A << 1, 1, 1e-20, 0;

		B << 0.5, 1;

		auto LU = A.lu();
		auto C = LU.solve(B);

		BOOST_CHECK(bertini::LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=bertini::MatrixSuccessCode::Success);

	}

	BOOST_AUTO_TEST_CASE(small_value_double)
	{
		BOOST_CHECK( bertini::IsSmallValue(std::complex<double>(1e-15,0)));
		BOOST_CHECK( bertini::IsSmallValue(std::complex<double>(1e-14,0)));

		BOOST_CHECK( bertini::IsSmallValue(std::complex<double>(-1e-15,0)));
		BOOST_CHECK( bertini::IsSmallValue(std::complex<double>(-1e-14,0)));

		BOOST_CHECK(!bertini::IsSmallValue(std::complex<double>(1e-10,0)));
		BOOST_CHECK(!bertini::IsSmallValue(std::complex<double>(1e2,0)));

		BOOST_CHECK(!bertini::IsSmallValue(std::complex<double>(-1e-10,0)));
		BOOST_CHECK(!bertini::IsSmallValue(std::complex<double>(-1e2,0)));
	}

	BOOST_AUTO_TEST_CASE(large_change_double)
	{
		BOOST_CHECK(bertini::IsLargeChange(1.0,1e-12));
		BOOST_CHECK(bertini::IsLargeChange(1e5,1e-7));
		BOOST_CHECK(!bertini::IsLargeChange(1e3,1e-4));
		BOOST_CHECK(!bertini::IsLargeChange(1e-15,1e-18));
	}


	BOOST_AUTO_TEST_CASE(small_value_multiprecision)
	{

		boost::multiprecision::mpfr_float::default_precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		boost::multiprecision::mpfr_float p = pow(mpfr_float(10),-mpfr_float(CLASS_TEST_MPFR_DEFAULT_DIGITS));

		BOOST_CHECK( bertini::IsSmallValue(bertini::complex(p,mpfr_float(0))));
		BOOST_CHECK( bertini::IsSmallValue(bertini::complex(-p,mpfr_float(0))));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(1e-15,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(1e-14,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(1e-10,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(1e2,0.0)));

		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(-1e-15,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(-1e-14,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(-1e-10,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex(-1e2,0.0)));
	}

	BOOST_AUTO_TEST_CASE(large_change_multiprecision)
	{
		boost::multiprecision::mpfr_float p = pow(mpfr_float(10),-mpfr_float(CLASS_TEST_MPFR_DEFAULT_DIGITS));

		BOOST_CHECK( bertini::IsLargeChange(mpfr_float(1.0),p));

		BOOST_CHECK( bertini::IsLargeChange(mpfr_float(1e16),mpfr_float(p*mpfr_float(1e16))));

		BOOST_CHECK( bertini::IsLargeChange(mpfr_float(-1e16),mpfr_float(p*mpfr_float(1e16))));
		BOOST_CHECK( bertini::IsLargeChange(mpfr_float(1e16),mpfr_float(-p*mpfr_float(1e16))));

		BOOST_CHECK(!bertini::IsLargeChange(mpfr_float(1.0),mpfr_float(1e-12)));
		BOOST_CHECK(!bertini::IsLargeChange(mpfr_float(-1.0),mpfr_float(1e-12)));
		BOOST_CHECK(!bertini::IsLargeChange(mpfr_float(1e5),mpfr_float(1e-7)));
		BOOST_CHECK(!bertini::IsLargeChange(mpfr_float(1e3),mpfr_float(1e-4)));
		BOOST_CHECK(!bertini::IsLargeChange(mpfr_float(1e-15),mpfr_float(1e-18)));
	}

	BOOST_AUTO_TEST_CASE(eigen_LU_partial_pivot_3x3)
	{
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(3,3);

		A << 0.000000010000000, 1.000000000000000,   1.000000000000000,
			 0 ,                1.000000000000000 ,  1.000000000000000,
			 1.000000000000000 ,1.000000000000000 ,  0;

		auto LU = A.lu();

		BOOST_CHECK(bertini::LUPartialPivotDecompositionSuccessful(LU.matrixLU())==bertini::MatrixSuccessCode::Success);
	}

BOOST_AUTO_TEST_SUITE_END()

	
	
	
	


