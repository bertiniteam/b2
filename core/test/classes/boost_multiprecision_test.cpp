#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpc.hpp>

using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>;
using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_on>;
using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_on>;
using mpc_complex = boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<0>, boost::multiprecision::et_on>;

namespace utf = boost::unit_test;

// BOOST_AUTO_TEST_SUITE(boost_multiprecision_tutorial)


// using mp_t = boost::multiprecision::debug_adaptor_t<boost::multiprecision::mpfr_float>;

// struct scoped_mpfr_precision
// {
// 	unsigned saved_digits10;
// 	scoped_mpfr_precision(unsigned digits10) : saved_digits10(mp_t::thread_default_precision())
// 	{
// 		mp_t::thread_default_precision(digits10);
// 	}

// 	~scoped_mpfr_precision()
// 	{
// 		mp_t::thread_default_precision(saved_digits10);
// 	}

// 	void reset(unsigned digits10)
// 	{
// 		mp_t::thread_default_precision(digits10);
// 	}

// 	void reset()
// 	{
// 		mp_t::thread_default_precision(saved_digits10);
// 	}
// };

// struct scoped_mpfr_precision_options
// {

// 	boost::multiprecision::variable_precision_options saved_options;

// 	scoped_mpfr_precision_options(boost::multiprecision::variable_precision_options opts) : saved_options(mp_t::thread_default_variable_precision_options())
// 	{
// 		mp_t::thread_default_variable_precision_options(opts);
// 	}

// 	~scoped_mpfr_precision_options()
// 	{
// 		mp_t::thread_default_variable_precision_options(saved_options);
// 	}

// 	void reset(boost::multiprecision::variable_precision_options opts)
// 	{
// 		mp_t::thread_default_variable_precision_options(opts);
// 	}
// };



// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_source)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_source_precision);


// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;
// 	std::cout << i << std::endl;
// 	mp_t f = i;

// 	BOOST_CHECK(f.precision()>mp_t::thread_default_precision());
// }


// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_target)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_target_precision);


// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;

// 	mp_t f(0);
// 	f = i;

// 	BOOST_CHECK_EQUAL(f.precision(),mp_t::default_precision());
// }


// BOOST_AUTO_TEST_CASE(precision_2_to_the_1000_minus_1_preserve_all)
// {
// 	scoped_mpfr_precision scope_1(mp_t::thread_default_precision());
// 	scoped_mpfr_precision_options scope_2(boost::multiprecision::variable_precision_options::preserve_all_precision);

// 	// calculate 2^1000 - 1:
// 	mpz_int i(1);
// 	i = i << 1000;
// 	i -= 1;

// 	mp_t f = i;

// 	BOOST_CHECK(f.precision()>mp_t::thread_default_precision());
// }

// BOOST_AUTO_TEST_SUITE_END() // boost_multiprecision_tutorial












BOOST_AUTO_TEST_SUITE(boost_multiprecision)

BOOST_AUTO_TEST_CASE(precision_complex_rational_add)
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a+b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_add_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b+a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_CASE(precision_complex_rational_sub, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a-b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_sub_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b-a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_CASE(precision_complex_rational_mul, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a*b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}

BOOST_AUTO_TEST_CASE(precision_complex_rational_mul_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b*a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_div, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = a/b;

	BOOST_CHECK_EQUAL(c.precision(),30);
}


BOOST_AUTO_TEST_CASE(precision_complex_rational_div_other_order, *utf::depends_on("boost_multiprecision/precision_complex_rational_add"))
{
	mpc_complex::default_precision(30);

	mpq_rational a(1,2);
	mpc_complex b(0,1);

	mpc_complex c = b/a;

	BOOST_CHECK_EQUAL(c.precision(),30);
}



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float_16digits)
{
	mpfr_float::default_precision(30);
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(b.precision(),16);
	BOOST_CHECK_EQUAL(result.precision(),16);
}

BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_complex_16digits, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_float_16digits"))
{
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_both_16digits, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_complex_16digits"))
{
	mpfr_float::default_precision(16);
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}

BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float20_complex16, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_both_16digits"))
{
	mpfr_float::default_precision(20);
	mpc_complex::default_precision(16);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),16);
}




BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float16_complex20, *utf::depends_on("boost_multiprecision/precision_complex_longlong_div_set_float20_complex16"))
{
	mpfr_float::default_precision(16);
	mpc_complex::default_precision(20);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),20);
}



BOOST_AUTO_TEST_CASE(precision_complex_longlong_div_set_float20_complex20)
{
	mpfr_float::default_precision(20);
	mpc_complex::default_precision(20);

	mpc_complex b(1,0);

	unsigned long long d{13};

	mpc_complex result{b/d};
	BOOST_CHECK_EQUAL(result.precision(),20);
}
































// shamelessly taken from the documentation for Boost.Multiprecision
// https://www.boost.org/doc/libs/1_82_0/libs/multiprecision/doc/html/boost_multiprecision/tut/variable.html


// one for all threads.
// as of 1.81, 1.82, there's a bug in Boost.Multiprecision, and the all threads version does NOT affect the current thread.
struct scoped_mpfr_precision_options_all_threads
{
   boost::multiprecision::variable_precision_options saved_options;
   scoped_mpfr_precision_options_all_threads(boost::multiprecision::variable_precision_options opts) : saved_options(mpfr_float::default_variable_precision_options())
   {
      mpfr_float::default_variable_precision_options(opts);
   }
   ~scoped_mpfr_precision_options_all_threads()
   {
      mpfr_float::default_variable_precision_options(saved_options);
   }
   void reset(boost::multiprecision::variable_precision_options opts)
   {
      mpfr_float::default_variable_precision_options(opts);
   }
};




// one for just this thread.
struct scoped_mpfr_precision_options_this_thread
{
   boost::multiprecision::variable_precision_options saved_options;
   scoped_mpfr_precision_options_this_thread(boost::multiprecision::variable_precision_options opts) : saved_options(mpfr_float::thread_default_variable_precision_options())
   {
      mpfr_float::thread_default_variable_precision_options(opts);
   }
   ~scoped_mpfr_precision_options_this_thread()
   {
      mpfr_float::thread_default_variable_precision_options(saved_options);
   }
   void reset(boost::multiprecision::variable_precision_options opts)
   {
      mpfr_float::thread_default_variable_precision_options(opts);
   }
};















BOOST_AUTO_TEST_CASE(precision_opts_all_threads_affect_current_thread)
{
	scoped_mpfr_precision_options_this_thread scoped_opts_this(boost::multiprecision::variable_precision_options::preserve_target_precision);

	scoped_mpfr_precision_options_all_threads scoped_opts_all(boost::multiprecision::variable_precision_options::preserve_source_precision);

	bool options_match = mpfr_float::default_variable_precision_options() == mpfr_float::thread_default_variable_precision_options();
	BOOST_CHECK(options_match && "there's a bug in Boost.Multiprecision in 1.82 and below: default_variable_precision_options() doesn't affect the current thread.  see https://github.com/boostorg/multiprecision/issues/551");
}










BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_all_threads_preserve_target_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread") )
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_target_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}





BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_all_threads_preserve_source_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_source_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}






BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_all_threads_preserve_component_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_component_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_all_threads_preserve_related_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_related_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_all_threads_preserve_all_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_all_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}





























BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_all_threads_preserve_target_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// // All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. The precision of other types (including related or component types - see preserve_component_precision/preserve_related_precision) contained within the expression are ignored. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_all_threads scoped_opts(boost::multiprecision::variable_precision_options::preserve_target_precision);
	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70); //
}
	




BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_all_threads_preserve_source_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");


	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);



	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and then rounded to the precision of the target variable upon assignment. The precision of other types (including related or component types - see preserve_component_precision/preserve_related_precision) contained within the expression are ignored. This option has the unfortunate side effect, that moves may become full deep copies. 

	scoped_mpfr_precision_options_all_threads   scoped_opts(boost::multiprecision::variable_precision_options::preserve_source_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_all_threads_preserve_component_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. If the expression contains component types then these are also considered when calculating the precision of the expression. Component types are the types which make up the two components of the number when dealing with interval or complex numbers. They are the same type as Num::value_type. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_all_threads   scoped_opts(boost::multiprecision::variable_precision_options::preserve_component_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}






BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_all_threads_preserve_related_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. If the expression contains component types then these are also considered when calculating the precision of the expression. In addition to component types, all related types are considered when evaluating the precision of the expression. Related types are considered to be instantiations of the same template, but with different parameters. So for example mpfr_float_100 would be a related type to mpfr_float, and all expressions containing an mpfr_float_100 variable would have at least 100 decimal digits of precision when evaluated as an mpfr_float expression. Moves, are true moves not copies.  

	scoped_mpfr_precision_options_all_threads   scoped_opts(boost::multiprecision::variable_precision_options::preserve_related_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}



BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_all_threads_preserve_all_precision, *utf::depends_on("boost_multiprecision/precision_opts_all_threads_affect_current_thread"))
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. In addition to component and related types, all types are considered when evaluating the precision of the expression. For example, if the expression contains an mpz_int, then the precision of the expression will be sufficient to store all of the digits in the integer unchanged. This option should generally be used with extreme caution, as it can easily cause unintentional precision inflation. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_all_threads   scoped_opts(boost::multiprecision::variable_precision_options::preserve_all_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}


































































BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_this_thread_preserve_target_precision)
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_target_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 30);
	
}





BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_this_thread_preserve_source_precision)
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_source_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}






BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_this_thread_preserve_component_precision)
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_component_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_this_thread_preserve_related_precision)
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_related_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_fresh_variable_this_thread_preserve_all_precision)
{



	mpfr_float::default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123");
	BOOST_CHECK_EQUAL(y.precision(), 60);



	mpfr_float::default_precision(30);
	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_all_precision);
	
	mpfr_float z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
	
}


















// this is expected to fail in Boost 1.81, 1.82, and presumably in earlier versions, due to a bug in Boost.Multiprecision in 1.81.  
// See https://github.com/boostorg/multiprecision/issues/551
// 
BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_this_thread_preserve_target_precision)
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// // All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. The precision of other types (including related or component types - see preserve_component_precision/preserve_related_precision) contained within the expression are ignored. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_this_thread scoped_opts(boost::multiprecision::variable_precision_options::preserve_target_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);
}
	




BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_this_thread_preserve_source_precision)
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");


	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);



	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and then rounded to the precision of the target variable upon assignment. The precision of other types (including related or component types - see preserve_component_precision/preserve_related_precision) contained within the expression are ignored. This option has the unfortunate side effect, that moves may become full deep copies. 

	scoped_mpfr_precision_options_this_thread   scoped_opts(boost::multiprecision::variable_precision_options::preserve_source_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}




BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_this_thread_preserve_component_precision)
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. If the expression contains component types then these are also considered when calculating the precision of the expression. Component types are the types which make up the two components of the number when dealing with interval or complex numbers. They are the same type as Num::value_type. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_this_thread   scoped_opts(boost::multiprecision::variable_precision_options::preserve_component_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}






BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_this_thread_preserve_related_precision)
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. If the expression contains component types then these are also considered when calculating the precision of the expression. In addition to component types, all related types are considered when evaluating the precision of the expression. Related types are considered to be instantiations of the same template, but with different parameters. So for example mpfr_float_100 would be a related type to mpfr_float, and all expressions containing an mpfr_float_100 variable would have at least 100 decimal digits of precision when evaluated as an mpfr_float expression. Moves, are true moves not copies.  

	scoped_mpfr_precision_options_this_thread   scoped_opts(boost::multiprecision::variable_precision_options::preserve_related_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}



BOOST_AUTO_TEST_CASE(arithmetic_precision_existing_variable_this_thread_preserve_all_precision)
{



	mpfr_float::thread_default_precision(50);
	mpfr_float x("0.01234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(x.precision(), 50);


	mpfr_float::thread_default_precision(60);
	mpfr_float y("0.12345678901234567890123456789012345678901234567890123456789");
	BOOST_CHECK_EQUAL(y.precision(), 60);


	// make a variable at precision 70
	mpfr_float::thread_default_precision(70);
	mpfr_float z("0");

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 70);


	mpfr_float::thread_default_precision(30);
	// then try to write into this existing variable with various policies


	// All expressions are evaluated at the precision of the highest precision variable within the expression, and that precision is preserved upon assignment. In addition to component and related types, all types are considered when evaluating the precision of the expression. For example, if the expression contains an mpz_int, then the precision of the expression will be sufficient to store all of the digits in the integer unchanged. This option should generally be used with extreme caution, as it can easily cause unintentional precision inflation. Moves, are true moves not copies. 

	scoped_mpfr_precision_options_this_thread   scoped_opts(boost::multiprecision::variable_precision_options::preserve_all_precision);

	
	z = x*y;

	BOOST_CHECK_EQUAL(x.precision(), 50);
	BOOST_CHECK_EQUAL(y.precision(), 60);
	BOOST_CHECK_EQUAL(z.precision(), 60);
}








































BOOST_AUTO_TEST_SUITE_END() // boost_multiprecision tests