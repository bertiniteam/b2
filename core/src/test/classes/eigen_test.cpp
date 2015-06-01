
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <boost/test/unit_test.hpp>


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>

#include "mpfr_complex.hpp"

#include <boost/timer/timer.hpp>





template <typename NumberType>
Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> kahan_matrix(unsigned int mat_size, NumberType c)
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






BOOST_AUTO_TEST_SUITE(kahan_matrix_solving_LU)



	BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_double) {
		unsigned int size = 10;
		srand(2);  rand();
		
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A = kahan_matrix(size, 0.285), B(size,size), C;
		
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
			kahan_matrix(size, boost::multiprecision::mpfr_float(0.285)), B(size,size), C;
		
		
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
		
		mpfr_matrix A = kahan_matrix(size, mpfr(0.285)), B(size,size), C;
		
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
		kahan_matrix(size, std::complex<double>(0.285)), B(size,size), C;
		
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
		kahan_matrix(size, bertini::complex("0.285","0.0")), B(size,size), C;
		
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
	
BOOST_AUTO_TEST_SUITE_END()

	
	
	
	


