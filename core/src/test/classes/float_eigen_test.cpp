


#include <boost/test/unit_test.hpp>

#include "float.hpp"


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>


BOOST_AUTO_TEST_SUITE(eigen_float_compatibility)





//	function [A] = kahan(c,n)
//	%
//	%       this function generates a scaled kahan matrix.
//	%
//	s = sqrt((1-c) * (1+c));
//	scale = 1;
//	A= -c * triu(ones(n)) + (1+c) * eye(n);
//	for j = 1:n
//	A(j,1:n) = A(j,1:n) * scale;
//	scale = scale * s;
//	end
//	%
//	%       scale the columns of kahan matrix.
//	%
//	for j = 1:n
//	A(:,j) = A(:,j)/j;
//	end
//	%
//	%       for n = 100 and c = 0.285, if one does the following in matlab:
//	%           A = kahan(c,n); [Q,S,V] = svd(A); [L,U,P] = lu(V);
//	%       then the matrices L and U are ill conditioned (~1e9). Hence
//	%       GEPP gives a large growth factor for orthogonal matrices as well.
//	%
//	%





BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_double_precision)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float c(0.285), s, scale(1);
	s = sqrt( (bertini::Float(1)-c) * (bertini::Float(1)+c) );
	
	int mat_size = 100;
	
	
	Eigen::Matrix<bertini::Float, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size), B(mat_size,mat_size), C;
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		for (int jj=ii; jj<mat_size; jj++)
		{
			A(ii,jj) = -c * bertini::Float(1);
		}
	}
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		A(ii,ii) += (bertini::Float(1)+c);
	}
	
	
	for (int jj=0; jj<mat_size; jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (int jj=0;jj<mat_size;jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(kk,jj)/=bertini::Float(jj+1);
		}
	}
	
	

	for (int ii=0; ii<mat_size; ii++){
		B(ii,ii) = -bertini::Float(1)/bertini::Float(ii+1);
	}
	
	C = A.lu().solve(B);
	
	//add statement on the value of C to actually test
}

	



BOOST_AUTO_TEST_CASE(solve_100x100_kahan_matrix_multiple_precision)
{
	boost::multiprecision::mpfr_float::default_precision(100);
	
	bertini::Float c(0.285), s, scale(1);
	s = sqrt(( bertini::Float(1)-c) * (bertini::Float(1)+c));
	
	int mat_size = 100;
	
	
	Eigen::Matrix<bertini::Float, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size), B(mat_size,mat_size), C;
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		for (int jj=ii; jj<mat_size; jj++)
		{
			A(ii,jj) = -c * bertini::Float(1);
		}
	}
	
	
	for (int ii=0; ii<mat_size; ii++)
	{
		A(ii,ii) += (bertini::Float(1)+c);
	}
	
	
	for (int jj=0; jj<mat_size; jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(jj,kk) *= scale;
		}
		scale *= s;
	}
	
	
	for (int jj=0;jj<mat_size;jj++){
		for (int kk=0;kk<mat_size;kk++){
			A(kk,jj)/=bertini::Float(jj+1);
		}
	}
	
	
	
	for (int ii=0; ii<mat_size; ii++){
		B(ii,ii) = -bertini::Float(1)/bertini::Float(ii+1);
	}
	
	C = A.lu().solve(B);
}





BOOST_AUTO_TEST_CASE(static_epsilon_function)
{
	
	boost::multiprecision::mpfr_float::default_precision(50);
	BOOST_CHECK_EQUAL(Eigen::NumTraits<bertini::Float>::epsilon(),bertini::Float(1e-49));
	
	
	boost::multiprecision::mpfr_float::default_precision(100);
	BOOST_CHECK_EQUAL(Eigen::NumTraits<bertini::Float>::epsilon(),bertini::Float(1e-99));
	
	
}




BOOST_AUTO_TEST_SUITE_END()


	


	
	
	
	
	
	
	
	


