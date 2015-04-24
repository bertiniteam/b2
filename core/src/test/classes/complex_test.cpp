
#include "complex.hpp"
#include <boost/test/unit_test.hpp>

#include <fstream>
#define private public

BOOST_AUTO_TEST_SUITE(complex_multiprecision_class)

BOOST_AUTO_TEST_CASE(complex_create_default_constructor)
{
	bertini::complex z;
}


BOOST_AUTO_TEST_CASE(complex_create_two_input_constructor)
{
	bertini::complex z("0.1","1.2");
}


BOOST_AUTO_TEST_CASE(complex_addition)
{
	bertini::complex z("0.1","1.2");
	bertini::complex v("0.2","1.3");
	
	bertini::complex r = z+v;
	BOOST_CHECK_EQUAL(r.real(),boost::multiprecision::mpfr_float("0.3"));
	BOOST_CHECK_EQUAL(r.imag(),boost::multiprecision::mpfr_float("2.5"));
}



BOOST_AUTO_TEST_CASE(complex_multiplication)
{
	boost::multiprecision::mpfr_float::default_precision(50);
	bertini::complex z("0.1000","1.2000");
	bertini::complex v("0.2000","1.3000");
	bertini::complex r = z*v;
	
	BOOST_CHECK(abs(r.real()-boost::multiprecision::mpfr_float("-1.54")) < 1e-50);
	BOOST_CHECK(abs(r.imag()-boost::multiprecision::mpfr_float("0.37"))< 1e-50);
	
}






BOOST_AUTO_TEST_CASE(complex_abs_of_i_is_one)
{
	bertini::complex i("0.0","1.0");
	boost::multiprecision::mpfr_float one = abs(i);
	std::cout << one << std::endl;
	BOOST_CHECK_EQUAL(one,boost::multiprecision::mpfr_float("1.0"));
}



BOOST_AUTO_TEST_CASE(complex_write_to_screen)
{
	bertini::complex i("0.0","1.0");
	
	std::cout << i.real() << " " << i.imag() << std::endl;
	std::cout << i << std::endl;
}








BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_with_comma)
{
	std::stringstream ss;
	ss << "(0.1234,-4.12397)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("0.1234"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("-4.12397"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_parens_no_comma)
{
	std::stringstream ss;
	ss << "(-3.651263418976498712e-2)";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("0"));
}


BOOST_AUTO_TEST_CASE(complex_get_from_stream_no_parens)
{
	std::stringstream ss;
	ss << "-3.651263418976498712e-2";
	
	bertini::complex z;
	ss >> z;
	
	BOOST_CHECK_EQUAL(z.real(),boost::multiprecision::mpfr_float("-3.651263418976498712e-2"));
	BOOST_CHECK_EQUAL(z.imag(),boost::multiprecision::mpfr_float("0"));
}



BOOST_AUTO_TEST_CASE(complex_eleventh_power)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	bertini::complex z("2.0","0.5"), w;
	
	w = pow(z,11);
	BOOST_CHECK_EQUAL(real(w), mpfr_float("-1319867") / mpfr_float("512") );
	BOOST_CHECK_EQUAL(imag(w),  mpfr_float("2529647") / mpfr_float("2048"));
	
	
}




BOOST_AUTO_TEST_CASE(complex_conjugation)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	mpfr_float::default_precision(50);
	
	bertini::complex z("2.0", sqrt(mpfr_float("2.0")));
	bertini::complex w = conj(z);
	
	BOOST_CHECK_EQUAL(real(w),mpfr_float("2.0"));
	BOOST_CHECK_EQUAL(imag(w),-sqrt(mpfr_float("2.0")));
}


BOOST_AUTO_TEST_CASE(complex_absolute_value)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	mpfr_float::default_precision(50);
	
	BOOST_CHECK_EQUAL(abs(bertini::complex("1.0","1.0")), sqrt(mpfr_float("2.0")));
	
}



BOOST_AUTO_TEST_CASE(complex_argument)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	mpfr_float::default_precision(50);
	
	bertini::complex z("1.0","1.0");
	
	BOOST_CHECK_EQUAL(arg(z), acos(mpfr_float("-1.0"))/mpfr_float("4.0"));
}








BOOST_AUTO_TEST_CASE(complex_serialization)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	mpfr_float::default_precision(50);
	
	
	bertini::complex z("1.23456", acos(mpfr_float("-1")));
	std::cout << z << std::endl;
	
	{
		std::ofstream fout("serialization_test_complex");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << z;
	}
	
	//	bertini::complex w;
	bertini::complex w;
	{
		std::ifstream fin("serialization_test_complex");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(real(z),real(w));
	BOOST_CHECK_EQUAL(imag(z),imag(w));
	
}





BOOST_AUTO_TEST_CASE(mpfr_float_serialization)
{
	using mpfr_float = boost::multiprecision::mpfr_float;
	mpfr_float::default_precision(50);
	
	
	mpfr_float q = acos( mpfr_float("-1.0") );
	
	
	{
		std::ofstream fout("serialization_test_mpfr");
	
		boost::archive::text_oarchive oa(fout);

		// write class instance to archive
		oa << q;
	}

	mpfr_float w;
	{
		std::ifstream fin("serialization_test_mpfr");
	
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> w;
	}
	
	BOOST_CHECK_EQUAL(q,w);
	
}


//// create and open a character archive for output
//std::ofstream ofs("float_archive_test");
//
//boost::multiprecision::mpfr_float::default_precision(16);
//
//bertini::Float a;
//a.MakeRandom();
//
//{
//	boost::archive::text_oarchive oa(ofs);
//	
//	// write class instance to archive
//	oa << a;
//	// archive and stream closed when destructors are called
//	}
//	
//	boost::multiprecision::mpfr_float::default_precision(100);
//	bertini::Float b;
//	{
//		// create and open an archive for input
//		std::ifstream ifs("float_archive_test");
//		boost::archive::text_iarchive ia(ifs);
//		// read class state from archive
//		ia >> b;
//		// archive and stream closed when destructors are called
//	}
//	
//	BOOST_CHECK_EQUAL(a,b);
//	BOOST_CHECK_EQUAL(b.Precision(),16);



BOOST_AUTO_TEST_SUITE_END()


