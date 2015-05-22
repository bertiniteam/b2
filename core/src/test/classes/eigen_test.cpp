
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <boost/test/unit_test.hpp>


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>



#include <boost/timer/timer.hpp>

BOOST_AUTO_TEST_SUITE(eigen_compatibility)











BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_pure_doubles)
{
	
	
	double c(0.285), s, scale(1);
	s = sqrt( (1-c) * (1+c) );
	
	int mat_size = 100;
	
	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size), B(mat_size,mat_size), C;
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		for (int jj=ii; jj<mat_size; jj++)
		{
			A(ii,jj) = -c * 1;
		}
	}
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		A(ii,ii) += (1+c);
	}
	
	
	for (int jj=0; jj<mat_size; jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (int jj=0;jj<mat_size;jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(kk,jj)/=jj+1;
		}
	}
	
	
	
	for (int ii=0; ii<mat_size; ii++){
		B(ii,ii) = -1/(ii+1);
	}
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "pure double time to solve:" << std::endl;
	//add statement on the value of C to actually test
}



BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_pure_mpfr_float_16)
{
	
	
	boost::multiprecision::mpfr_float::default_precision(16);
	
	boost::multiprecision::mpfr_float c(0.285), s, scale(1);
	s = sqrt(( boost::multiprecision::mpfr_float(1)-c) * (boost::multiprecision::mpfr_float(1)+c));
	
	int mat_size = 100;
	
	
	Eigen::Matrix<boost::multiprecision::mpfr_float, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size), B(mat_size,mat_size), C;
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		for (int jj=ii; jj<mat_size; jj++)
		{
			A(ii,jj) = -c * boost::multiprecision::mpfr_float(1);
		}
	}
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		A(ii,ii) += (boost::multiprecision::mpfr_float(1)+c);
	}
	
	
	for (int jj=0; jj<mat_size; jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (int jj=0;jj<mat_size;jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(kk,jj)/=boost::multiprecision::mpfr_float(jj+1);
		}
	}
	
	
	
	for (int ii=0; ii<mat_size; ii++){
		B(ii,ii) = -boost::multiprecision::mpfr_float(1)/boost::multiprecision::mpfr_float(ii+1);
	}
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "mpfr_float pure 16 time to solve:" << std::endl;
}




BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_pure_mpfr_float_100)
{
	
	
	boost::multiprecision::mpfr_float::default_precision(100);
	
	boost::multiprecision::mpfr_float c(0.285), s, scale(1);
	s = sqrt(( boost::multiprecision::mpfr_float(1)-c) * (boost::multiprecision::mpfr_float(1)+c));
	
	int mat_size = 100;
	
	
	Eigen::Matrix<boost::multiprecision::mpfr_float, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size), B(mat_size,mat_size), C;
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		for (int jj=ii; jj<mat_size; jj++)
		{
			A(ii,jj) = -c * boost::multiprecision::mpfr_float(1);
		}
	}
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		A(ii,ii) += (boost::multiprecision::mpfr_float(1)+c);
	}
	
	
	for (int jj=0; jj<mat_size; jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (int jj=0;jj<mat_size;jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(kk,jj)/=boost::multiprecision::mpfr_float(jj+1);
		}
	}
	
	
	
	for (int ii=0; ii<mat_size; ii++){
		B(ii,ii) = -boost::multiprecision::mpfr_float(1)/boost::multiprecision::mpfr_float(ii+1);
	}
	
	boost::timer::auto_cpu_timer t;
	C = A.lu().solve(B);
	
	std::cout << "mpfr_float pure 100 time to solve:" << std::endl;
}



	
	
	
BOOST_AUTO_TEST_SUITE_END()


	


	
	
	
	
	
	
	
	


