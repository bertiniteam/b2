
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <boost/test/unit_test.hpp>


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>

#include "complex.hpp"

#include <boost/timer/timer.hpp>






template <typename numbertype>
Eigen::Matrix<numbertype, Eigen::Dynamic, Eigen::Dynamic> kahan_matrix(unsigned int mat_size, numbertype c)
{
	numbertype s, scale(1.0);
	s = sqrt( (numbertype(1.0)-c) * (numbertype(1.0)+c) );
	
	
	Eigen::Matrix<numbertype, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size);
	
	
	for (unsigned int ii=0; ii<mat_size; ii++) {
		for (unsigned int jj=ii; jj<mat_size; jj++) {
			A(ii,jj) = -c * numbertype(1.0);
		}
	}
	
	
	for (unsigned int ii=0; ii<mat_size; ii++) {
		A(ii,ii) += numbertype(1)+c;
	}
	
	
	for (unsigned int jj=0; jj<mat_size; jj++){
		for (unsigned int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (unsigned int jj=0;jj<mat_size;jj++){
		for (unsigned int kk=0;kk<mat_size;kk++){
			A(kk,jj)/= numbertype(jj) + numbertype(1);
		}
	}
	
	return A;
}






BOOST_AUTO_TEST_SUITE(kahan_matrix_solving_LU)



BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_double) {
	
	srand(2);  rand();
	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A = kahan_matrix(100, 0.285), B(100,100), C;
	
	for (int ii=0; ii<100; ii++)
		B(ii,ii) = -1/(ii+1) + double(rand()) /  RAND_MAX;
	
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "pure double time to solve:" << std::endl;
	//add statement on the value of C to actually test

}



BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_16)
{
	
	
	boost::multiprecision::mpfr_float::default_precision(16);
	
	srand(2);  rand();
	
	Eigen::Matrix<boost::multiprecision::mpfr_float, Eigen::Dynamic, Eigen::Dynamic> A =
		kahan_matrix(100, boost::multiprecision::mpfr_float(0.285)), B(100,100), C;
	
	for (int ii=0; ii<100; ii++)
		B(ii,ii) = -1/(ii+1) + boost::multiprecision::mpfr_float(rand()) /  boost::multiprecision::mpfr_float(RAND_MAX);
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "mpfr_float pure 16 time to solve:" << std::endl;
}




BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_mpfr_float_100)
{
	
	
	boost::multiprecision::mpfr_float::default_precision(100);
	
	srand(2);  rand();
	
	Eigen::Matrix<boost::multiprecision::mpfr_float, Eigen::Dynamic, Eigen::Dynamic> A =
	kahan_matrix(100, boost::multiprecision::mpfr_float(0.285)), B(100,100), C;
	
	for (int ii=0; ii<100; ii++)
		B(ii,ii) = -1/(ii+1) + boost::multiprecision::mpfr_float(rand()) /  boost::multiprecision::mpfr_float(RAND_MAX);
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "mpfr_float pure 100 time to solve:" << std::endl;
	
}


	
	
	
//BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_bertinicomplex_100)
//{
//	
//	
//	boost::multiprecision::mpfr_float::default_precision(100);
//	
//	srand(2);  rand();
//	
//	Eigen::Matrix<bertini::complex, Eigen::Dynamic, Eigen::Dynamic> A =
//	kahan_matrix(100, bertini::complex("0.285","0.0")), B(100,100), C;
//	
//	for (int ii=0; ii<100; ii++)
//		B(ii,ii) = bertini::complex(-1/(ii+1) + double(rand()) /  double(RAND_MAX), 0.0);
//	
//	boost::timer::auto_cpu_timer t;
//	C = A.lu().solve(B);
//	
//	std::cout << "bertini::complex precision 100 time to solve:" << std::endl;
//}

	
BOOST_AUTO_TEST_SUITE_END()

	
	
	
	


