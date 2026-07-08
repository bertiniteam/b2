//This file is part of Bertini 2.
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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire




#include <boost/test/unit_test.hpp>





#include "bertini2/eigen_extensions.hpp"

#include <Eigen/Dense>
#include <Eigen/LU>

#include "externs.hpp"




BOOST_AUTO_TEST_SUITE(num_traits_checking)

using real_mp = bertini::real_mp;

// this test assures that the Eigen::NumTraits defined in eigen_extensions.hpp is correctly found during template instantiation, and that the Real type it defines is actually real_mp.  If it were not, then we would be unable to make an expression of ::Real type in variable q, because it would be an expression, not a populatable number of type real_mp.
BOOST_AUTO_TEST_CASE(expressions_of_mpfr_floats)
{
	using NumT = bertini::real_mp;

	NumT a {1}, b{2}, c{3};

	Eigen::NumTraits<decltype(a*a + b*b/ c)>::Real q{0};

}

BOOST_AUTO_TEST_CASE(size_object_sensible_vec)
{
	bertini::Vec<bertini::complex_dbl> v(3);
	BOOST_CHECK_EQUAL(v.rows(),3);
	BOOST_CHECK_EQUAL(v.cols(),1);

	BOOST_CHECK(!bertini::IsEmpty(v));
}


BOOST_AUTO_TEST_CASE(size_object_sensible_mat)
{
	bertini::Mat<bertini::complex_dbl> v(3,4);
	BOOST_CHECK_EQUAL(v.rows(),3);
	BOOST_CHECK_EQUAL(v.cols(),4);

	BOOST_CHECK(!bertini::IsEmpty(v));
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(kahan_matrix_solving_LU)

using real_mp = bertini::real_mp;
using bertini::KahanMatrix;

	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_double) {
		unsigned int size = 10;
		srand(2);  rand();
		
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A = KahanMatrix(size, 0.285), B(size,size), C;
		
		for (unsigned int ii=0; ii<size; ii++)
			for (unsigned int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -1.0/(ii+1) + double(rand()) /  RAND_MAX : B(ii,jj) = 0;
		
		C = A.lu().solve(B);
		//add statement on the value of C to actually test

	}



	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_16)
	{
		
		unsigned int size = 10;
		bertini::DefaultPrecision(16);
		
		srand(2);  rand();
		
		Eigen::Matrix<bertini::real_mp, Eigen::Dynamic, Eigen::Dynamic> A =
			KahanMatrix(size, bertini::real_mp(0.285)), B(size,size), C;
		
		
		for (unsigned int ii=0; ii<size; ii++)
			for (unsigned int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -bertini::real_mp(1)/(ii+1) + bertini::real_mp(rand()) /  bertini::real_mp(RAND_MAX) : B(ii,jj) = 0;

		C = A.lu().solve(B);
	}




	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_100)
	{
		unsigned int size = 10;
		srand(2);  rand();
		
		using mpfr = bertini::real_mp;
		using mpfr_matrix = Eigen::Matrix<mpfr, Eigen::Dynamic, Eigen::Dynamic>;

		bertini::DefaultPrecision(100);
		
		mpfr_matrix A = KahanMatrix(size, mpfr(0.285)), B(size,size), C;
		
		for (unsigned int ii=0; ii<size; ii++){
			for (unsigned int jj=0; jj<size; jj++){
				(jj!=ii) ? B(ii,jj) = -mpfr(1)/(ii+1) + mpfr(rand()) /  mpfr(RAND_MAX) : B(ii,jj) = mpfr(0.0);
			}
		}
		C = A.lu().solve(B);

	}


	
	

	
	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_standardcomplex)
	{
		unsigned int size = 10;
		srand(2);  rand();
		
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A =
		KahanMatrix(size, std::complex<double>(0.285)), B(size,size), C;
		
		for (unsigned int ii=0; ii<size; ii++)
			for (unsigned int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = -1.0/(ii+1) + double(rand()) / double(RAND_MAX) : B(ii,jj) = 0;
		
		
		C = A.lu().solve(B);
	}
	
	
	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_bertinicomplex_100)
	{
		
		unsigned int size = 10;
		bertini::DefaultPrecision(100);
		
		srand(2);  rand();
		
		Eigen::Matrix<bertini::complex_mp, Eigen::Dynamic, Eigen::Dynamic> A =
		KahanMatrix(size, bertini::complex_mp("0.285","0.0")), B(size,size), C;
		
		for (unsigned int ii=0; ii<size; ii++)
			for (unsigned int jj=0; jj<size; jj++)
				jj!=ii? B(ii,jj) = bertini::complex_mp( bertini::complex_mp(-1)/bertini::complex_mp(ii+1) + bertini::complex_mp(rand()) / bertini::complex_mp(RAND_MAX)) : B(ii,jj) = bertini::complex_mp(0);
		
		C = A.lu().solve(B);
	}

		
		
		
		
//		BOOST_AUTO_TEST_CASE(mpfr_float_num_traits){
//			
//			std::cout << Eigen::NumTraits<bertini::real_mp>::highest() << std::endl;
//			std::cout << Eigen::NumTraits<bertini::real_mp>::lowest() << std::endl;
//			std::cout << Eigen::NumTraits<bertini::real_mp>::dummy_precision() << std::endl;
//			std::cout << Eigen::NumTraits<bertini::real_mp>::epsilon() << std::endl;
//			
//		}
	

	BOOST_AUTO_TEST_CASE(eigen_partial_pivot_solve_singular_matrix)
	{

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(2,2), B(2,1);

		A << 1, 1, 0, 0;

		B << 0.5, 1;

		auto LU = A.lu();


		[[maybe_unused]] auto C = LU.solve(B);

	}


	BOOST_AUTO_TEST_CASE(eigen_partial_pivot_solve_near_singular_matrix_double)
	{

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(2,2), B(2,1);

		A << 1, 1, 1e-20, 0;

		B << 0.5, 1;

		auto LU = A.lu();
		[[maybe_unused]] auto C = LU.solve(B);

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

		bertini::DefaultPrecision(30);

		bertini::real_mp p = pow(real_mp(10),-real_mp(30));

		BOOST_CHECK( bertini::IsSmallValue(bertini::complex_mp(p,real_mp(0))));
		BOOST_CHECK( bertini::IsSmallValue(bertini::complex_mp(-p,real_mp(0))));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(1e-15,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(1e-14,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(1e-10,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(1e2,0.0)));

		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(-1e-15,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(-1e-14,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(-1e-10,0.0)));
		BOOST_CHECK(!bertini::IsSmallValue(bertini::complex_mp(-1e2,0.0)));
	}

	BOOST_AUTO_TEST_CASE(large_change_multiprecision)
	{
		bertini::DefaultPrecision(30);

		bertini::real_mp p = pow(real_mp(10),-real_mp(30));

		BOOST_CHECK( bertini::IsLargeChange(real_mp(1.0),p));

		BOOST_CHECK( bertini::IsLargeChange(real_mp(1e16),real_mp(p*real_mp(1e16))));

		BOOST_CHECK( bertini::IsLargeChange(real_mp(-1e16),real_mp(p*real_mp(1e16))));
		BOOST_CHECK( bertini::IsLargeChange(real_mp(1e16),real_mp(-p*real_mp(1e16))));

		BOOST_CHECK(!bertini::IsLargeChange(real_mp(1.0),real_mp(1e-12)));
		BOOST_CHECK(!bertini::IsLargeChange(real_mp(-1.0),real_mp(1e-12)));
		BOOST_CHECK(!bertini::IsLargeChange(real_mp(1e5),real_mp(1e-7)));
		BOOST_CHECK(!bertini::IsLargeChange(real_mp(1e3),real_mp(1e-4)));
		BOOST_CHECK(!bertini::IsLargeChange(real_mp(1e-15),real_mp(1e-18)));
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


	BOOST_AUTO_TEST_CASE(eigen_norm_of_vector)
	{
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(1,3);
		A << 1, 2, 3;
		[[maybe_unused]] double n = A.norm();
	}

	BOOST_AUTO_TEST_CASE(dot_product_with_mpfr_type)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, 3, 1> v(data_type(2),data_type(4),data_type(3));
		Eigen::Matrix<data_type, 3, 1> w(data_type(1),data_type(2),data_type(-1));

		data_type result = v.dot(w);
		data_type exact(7);
		
		BOOST_CHECK_EQUAL(result, exact);
		
	}


	BOOST_AUTO_TEST_CASE(svd_with_mpfr_type)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		// this breaks with boost multiprecision et_on with eigen 3.2.7.
		Eigen::JacobiSVD<Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		
		
	}


	// scalar multiplication with various types

	BOOST_AUTO_TEST_CASE(scalar_multiplication_mpfr_mpfr)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		data_type a(1);
		// this breaks with boost multiprecision et_on with eigen 3.2.7.
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B = a*A;
		B = A*a;
	}


	BOOST_AUTO_TEST_CASE(scalar_multiplication_mpfr_int)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		int a(1);

		[[maybe_unused]] auto b = a*q;

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B = a*A;
		B = A*a;
	}


	BOOST_AUTO_TEST_CASE(scalar_multiplication_mpfr_long)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		long a(1);

		[[maybe_unused]] auto b = a*q;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B = a*A;
		B = A*a;
	}

	BOOST_AUTO_TEST_CASE(scalar_multiplication_mpfr_mpz_int)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		bertini::mpz_int a(1);

		[[maybe_unused]] auto b = a*q;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B = a*A;
		B = A*a;
	}



	// self multiplication


	BOOST_AUTO_TEST_CASE(self_multiplication_dbl_int)
	{
		
		using data_type = bertini::complex_dbl;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		int a(1);

		A*=a;
	}


	BOOST_AUTO_TEST_CASE(self_multiplication_mpfr_mpfr)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		data_type a(1);

		A*=a;
	}


	BOOST_AUTO_TEST_CASE(self_multiplication_mpfr_int)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		int a(1);

		[[maybe_unused]] auto b = a*q;

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		A*=a;
	}


	BOOST_AUTO_TEST_CASE(self_multiplication_mpfr_long)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		long a(1);

		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		A*=a;
	}

	BOOST_AUTO_TEST_CASE(self_multiplication_mpfr_mpz_int)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using data_type = bertini::complex_mp;
		
		data_type q(1);
		bertini::mpz_int a(1);
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		
		A*=a;
	}

	BOOST_AUTO_TEST_CASE(change_precision_mpfr_float)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using bertini::Precision;
		using data_type = bertini::real_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);

		Precision(A,100);
		BOOST_CHECK_EQUAL(A(0,0).precision(),100);
	}

	BOOST_AUTO_TEST_CASE(change_precision_mpfr_complex)
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		using bertini::Precision;
		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);

		Precision(A,100);
		BOOST_CHECK_EQUAL(A(0,0).precision(),100);
	}



	BOOST_AUTO_TEST_CASE(change_precision_mpfr_complex2)
	{
		bertini::DefaultPrecision(50);

		using bertini::Precision;
		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(2,2);
		A << data_type(2), data_type(1), data_type(1), data_type(2);

		Precision(A,100);
		BOOST_CHECK_EQUAL(A(0,0).precision(),100);

		auto new_prec = 50-10;

		bertini::DefaultPrecision(static_cast<unsigned int>(new_prec));
		Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B(2,2);
		
		B = A; // assignment preserves precision of source
		BOOST_CHECK((A-B).norm() < 1e-38);
		BOOST_CHECK_EQUAL(Precision(B),100);

	} 

	BOOST_AUTO_TEST_CASE(copy_matrix)
	{
		bertini::DefaultPrecision(50);

		using bertini::Precision;
		using data_type = bertini::complex_mp;
		
		Eigen::Matrix<data_type, Eigen::Dynamic, 1> A(4);
		A << data_type(2), data_type(1), data_type(1), data_type(2);
		Precision(A,100);
		BOOST_CHECK_EQUAL(A(0).precision(),100);

		auto new_prec = 50-10;

		bertini::DefaultPrecision(static_cast<unsigned int>(new_prec));
		Eigen::Matrix<data_type, Eigen::Dynamic, 1> B(4);

		B = A;
		BOOST_CHECK((A-B).norm() < 1e-38);
		BOOST_CHECK_EQUAL(Precision(B),100);

	}

BOOST_AUTO_TEST_SUITE_END()


// issue #304: IsDistinct / IsSamePoint -- the infinity-norm tolerance point comparison.
BOOST_AUTO_TEST_SUITE(point_distinctness)

BOOST_AUTO_TEST_CASE(is_distinct_infinity_norm)
{
	using namespace bertini;

	Vec<complex_mp> a(2), b(2), c(2);
	a << complex_mp(1), complex_mp(2);
	b << complex_mp(1), complex_mp(2) + complex_mp("1e-9");   // 1e-9 away in one coordinate
	c << complex_mp(1), complex_mp(3);                        // 1 away

	BOOST_CHECK(!IsDistinct(a, b, 1e-6));    // within tol -> same point
	BOOST_CHECK( IsSamePoint(a, b, 1e-6));
	BOOST_CHECK( IsDistinct(a, c, 1e-6));    // far apart -> distinct
	BOOST_CHECK( IsDistinct(a, b, 1e-12));   // tighter tol -> distinct

	Vec<complex_mp> shorter(1);
	shorter << complex_mp(1);
	BOOST_CHECK(IsDistinct(a, shorter, 1e-6));   // different length -> distinct, never same

	// the same predicate works on plain doubles (so callers need not think about the number type)
	Vec<double> da(2), db(2);
	da << 1.0, 2.0;
	db << 1.0, 2.0 + 1e-9;
	BOOST_CHECK(!IsDistinct(da, db, 1e-6));
	BOOST_CHECK( IsDistinct(da, db, 1e-12));
}

BOOST_AUTO_TEST_SUITE_END()

	
	
	
	


