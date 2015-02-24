
#include <boost/test/unit_test.hpp>





#define private public
#include "float.hpp"













// start of the test suite using bertini::Float
BOOST_FIXTURE_TEST_SUITE( float_basics, bertini::Float )



/**
 an empty constructed bertini::Float should have 0 precision, and two nullptr unique_ptr's.
 */
BOOST_AUTO_TEST_CASE(empty_has_both_fields_null)
{
	BOOST_CHECK_EQUAL(current_precision_,boost::multiprecision::mpfr_float::default_precision());
	
	if (boost::multiprecision::mpfr_float::default_precision()<17) {
		BOOST_CHECK(number_as_double_!=nullptr);
		BOOST_CHECK(number_as_multiple_precision_==nullptr);
		BOOST_CHECK_EQUAL(*number_as_double_,0);
	}
	else{
		BOOST_CHECK(number_as_double_==nullptr);
		BOOST_CHECK(number_as_multiple_precision_!=nullptr);
		BOOST_CHECK_EQUAL(*number_as_multiple_precision_,boost::multiprecision::mpfr_float("0.0"));
	}
	
}



BOOST_AUTO_TEST_CASE(change_precision_to_double_has_null_mp)
{
	ChangePrecision(16);
	
	BOOST_CHECK_EQUAL(current_precision_,16);
	BOOST_CHECK(number_as_double_!=nullptr);
	BOOST_CHECK(number_as_multiple_precision_==nullptr);
	
}


BOOST_AUTO_TEST_CASE(change_precision_to_mp_has_null_double50)
{
	ChangePrecision(50);
	
	BOOST_CHECK_EQUAL(current_precision_,50);
	BOOST_CHECK(number_as_double_==nullptr);
	BOOST_CHECK(number_as_multiple_precision_!=nullptr);
	BOOST_CHECK_EQUAL(number_as_multiple_precision_->precision(),50);
	
	
}



BOOST_AUTO_TEST_CASE(change_precision_to_mp_has_null_double100)
{
	ChangePrecision(100);
	
	BOOST_CHECK_EQUAL(current_precision_,100);
	BOOST_CHECK(number_as_double_==nullptr);
	BOOST_CHECK(number_as_multiple_precision_!=nullptr);
	BOOST_CHECK_EQUAL(number_as_multiple_precision_->precision(),100);
}



BOOST_AUTO_TEST_CASE(make_random_multiple_precision_correct_precision30)
{
	ChangePrecision(30);
	MakeRandom();
	BOOST_CHECK_EQUAL(number_as_multiple_precision_->precision(),30);
}




BOOST_AUTO_TEST_CASE(make_random_multiple_precision_correct_precision50)
{
	ChangePrecision(50);
	MakeRandom();
	BOOST_CHECK_EQUAL(number_as_multiple_precision_->precision(),50);
}

BOOST_AUTO_TEST_CASE(make_random_multiple_precision_correct_precision5000)
{
	ChangePrecision(5000);
	MakeRandom();
	BOOST_CHECK(number_as_multiple_precision_->precision()>=5000);
}

BOOST_AUTO_TEST_SUITE_END()
// the end of the suite using a bertini::Float as a fixture










/*******************
 **
 **  exercise constructors
 **
 *********************/



BOOST_AUTO_TEST_SUITE(Float)


BOOST_AUTO_TEST_CASE(Float_constructor_from_number)
{
	
	boost::multiprecision::mpfr_float::default_precision(50);
	
	bertini::Float a("2.0",30);
	
	BOOST_CHECK_EQUAL(a.Precision(),30);
	
	bertini::Float b(a);
	
	BOOST_CHECK_EQUAL(b.Precision(),30);
	BOOST_CHECK_EQUAL(*a.number_as_multiple_precision_,*b.number_as_multiple_precision_);
	
}

/**
 changing the default precision in boost::multiprecision::mpfr_float actually changes that of a default constructed bertini::Float
 */
BOOST_AUTO_TEST_CASE(change_default_precision)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	boost::multiprecision::mpfr_float a;
	BOOST_CHECK_EQUAL(a.precision(),16);
	
	bertini::Float a_bertini;
	
	
	boost::multiprecision::mpfr_float::default_precision(75);
	boost::multiprecision::mpfr_float b;
	BOOST_CHECK_EQUAL(b.precision(),75);
	
	
	boost::multiprecision::mpfr_float::default_precision(4962);
	boost::multiprecision::mpfr_float c;
	BOOST_CHECK_EQUAL(c.precision(),4962);
	
}




/**
 create a bertini::Float from a string, but the default is higher, and despite the fact you put in a short number, you get a higher-precision number, because that's the default.
 */
BOOST_AUTO_TEST_CASE(Float_create_from_string_multiple_precision0)
{
	boost::multiprecision::mpfr_float::default_precision(1000);
	bertini::Float a("0.1234");
	BOOST_CHECK_EQUAL(a.Precision(),1000);
}


/**
 make a Float which is longer than default precision, yet it is constructed as the default precision.
 */
BOOST_AUTO_TEST_CASE(Float_create_from_string_multiple_precision1)
{
	boost::multiprecision::mpfr_float::default_precision(1000);
	
	bertini::Float a("0.12345678901234567890123456789");
	BOOST_CHECK_EQUAL(a.Precision(),1000);
}

/**
 precision of the desired string and the default match
 */
BOOST_AUTO_TEST_CASE(Float_create_from_string_multiple_precision2)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	bertini::Float a("-0.10000000000000000000000000000");
	BOOST_CHECK_EQUAL(a.Precision(),30);
}


/**
 
 */
BOOST_AUTO_TEST_CASE(Float_create_from_string_multiple_precision3)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	bertini::Float a("1e30");
	BOOST_CHECK_EQUAL(a.Precision(),30);
}


/**
 exercise the (string,digits10) constructor, for constructing a bertini::Float with precision different from the current default.
 */
BOOST_AUTO_TEST_CASE(Float_create_from_string_precision_different_from_default)
{
	boost::multiprecision::mpfr_float::default_precision(50);
	
	bertini::Float a("-0.10000000000000000000000000000",30);
	BOOST_CHECK_EQUAL(a.number_as_multiple_precision_->precision(),30);
	
}




/**
 round-trip (or convert-and-come-back) a bertini::Float into a string and back again, from one which was in double form.
 */
BOOST_AUTO_TEST_CASE(Float_round_tripping_double)
{
	srand(time(NULL));
	double burnone = double(rand());
	
	bertini::Float a;
	a.ChangePrecision(16);
	a.MakeRandom();
	
	std::stringstream converter;
	converter.precision(a.Precision());
	converter << a;
	bertini::Float b;  b.ChangePrecision(16);
	converter >> b;
	
	
	BOOST_CHECK_EQUAL(a,b);
}

/**
 round-trip (or convert-and-come-back) a bertini::Float into a string and back again, from one which was in multiple precision form.
 
 note that you have to have one more digit of precision on the stream for conversion, than in the bertini::Float itself.  I find this weird.  
 
 TODO: find out why this is.
 */
BOOST_AUTO_TEST_CASE(Float_round_tripping_multiple_precision)
{
	
	bertini::Float a;
	a.ChangePrecision(50);
	a.MakeRandom();
	
	std::stringstream converter;
	converter.precision(a.Precision()+1);
	converter << a;
	bertini::Float b;
	b.ChangePrecision(50);
	converter >> b;
	

	BOOST_CHECK_EQUAL(a,b);
}


/**
 find the absolute value of a bertini::Float which was in double form.
 */
BOOST_AUTO_TEST_CASE(Float_absolute_value_double)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float a(-0.1);
	
	BOOST_CHECK_EQUAL( abs(a) , bertini::Float(0.1));
	
}


/**
 check the absolute value of a bertini::Float which was in multiple precision
 */
BOOST_AUTO_TEST_CASE(Float_absolute_value_multiple_precision)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a("-0.1");
	BOOST_CHECK_EQUAL( abs(a) , bertini::Float(0.1));
	
}









BOOST_AUTO_TEST_CASE(Float_add_both_double_check_correct_precision_and_value)
{
	
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float a(1.0), b(2.0);
	
	bertini::Float c = a+b;
	
	BOOST_CHECK_EQUAL(c.Precision(),16);
	BOOST_CHECK_EQUAL(*c.number_as_double_,3.0);
	
}



BOOST_AUTO_TEST_CASE(Float_add_both_multiple_precision_check_correct_precision_and_value)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a("1.0",30), b("2.0");
	boost::multiprecision::mpfr_float expected_result("3.0");
	bertini::Float c = a+b;
	
	BOOST_CHECK_EQUAL(c.Precision(),30);
	BOOST_CHECK_EQUAL(*c.number_as_multiple_precision_,expected_result);
}



BOOST_AUTO_TEST_CASE(Float_add_mixed_double_and_multiple_precision_throws)
{
	boost::multiprecision::mpfr_float::default_precision(50);
	bertini::Float a(0), b("1.2"), c;
	
	BOOST_CHECK_NO_THROW(c = a+b);
	BOOST_CHECK_EQUAL(*c.number_as_multiple_precision_,boost::multiprecision::mpfr_float("1.2"));
}


BOOST_AUTO_TEST_CASE(Float_add_differing_multiple_precision_check_throws)
{
	bertini::Float a, b, c;
	a.ChangePrecision(30);
	b.ChangePrecision(70);
	
	a.MakeRandom();
	b.MakeRandom();
	
	BOOST_CHECK_THROW(c = a+b,std::exception);
}


BOOST_AUTO_TEST_CASE(Float_add_same_multiple_precision_check_no_throw)
{
	bertini::Float a, b, c;
	
	a.ChangePrecision(50);
	b.ChangePrecision(50);
	
	a.MakeRandom();
	b.MakeRandom();
	
	BOOST_CHECK_NO_THROW(c = a+b);
}














BOOST_AUTO_TEST_CASE(Float_subtract_both_double_check_correct_precision_and_value)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float a(1.0), b(2.0);
	
	bertini::Float c = a-b;
	
	BOOST_CHECK_EQUAL(c.Precision(),16);
	BOOST_CHECK_EQUAL( (*c.number_as_double_),-1.0);
	
}


BOOST_AUTO_TEST_CASE(Float_subtract_both_double_check_correct_precision_and_value_2)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a(1.0), b(2.0);
	
	bertini::Float c = a-b;
	
	BOOST_CHECK_EQUAL(c.Precision(),30);
	BOOST_CHECK_EQUAL( *c.number_as_multiple_precision_,boost::multiprecision::mpfr_float("-1.0"));
	
}



BOOST_AUTO_TEST_CASE(Float_subtract_both_multiple_precision_check_correct_precision_and_value)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a("1.0",30), b("2.0");
	boost::multiprecision::mpfr_float expected_result("-1.0");
	bertini::Float c = a-b;
	
	BOOST_CHECK_EQUAL(c.Precision(),30);
	BOOST_CHECK_EQUAL(*c.number_as_multiple_precision_,expected_result);
}



BOOST_AUTO_TEST_CASE(Float_subtract_differing_multiple_precision_throws)
{
	bertini::Float a("1.323",30), b("2.43455",70), c;

	
	BOOST_CHECK_THROW(c = a-b,std::exception);
}


BOOST_AUTO_TEST_CASE(Float_subtract_same_multiple_precision_no_throws)
{
	bertini::Float a("1.323",50), b("2.43455",50), c;
	
	
	BOOST_CHECK_NO_THROW(c = a-b);
}


















BOOST_AUTO_TEST_CASE(Float_multiply_both_double_check_correct_precision)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	bertini::Float a(1.5), b(2.0);
	
	bertini::Float c = a*b;
	
	BOOST_CHECK_EQUAL(c.Precision(),16);
	BOOST_CHECK_EQUAL(*c.number_as_double_,3.0);
	
}



BOOST_AUTO_TEST_CASE(Float_multiply_both_multiple_precision_check_correct_precision_and_value)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a("1.5",30), b("2.0");
	boost::multiprecision::mpfr_float expected_result("3.0");
	bertini::Float c = a*b;
	
	BOOST_CHECK_EQUAL(c.Precision(),30);
	BOOST_CHECK_EQUAL(*c.number_as_multiple_precision_,expected_result);
}




BOOST_AUTO_TEST_CASE(Float_multiply_both_multiple_precision_check_mixed_throws)
{
	bertini::Float a("1.323",30), b("2.43455",70), c;
	
	BOOST_CHECK_THROW(c = a*b,std::exception);
}








BOOST_AUTO_TEST_CASE(Float_divide_both_double_check_correct_precision)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	bertini::Float a(1.5), b(2.0);
	
	bertini::Float c = a/b;
	
	BOOST_CHECK_EQUAL(c.Precision(),16);
	BOOST_CHECK_EQUAL(*c.number_as_double_,0.75);
	
}


BOOST_AUTO_TEST_CASE(Float_divide_both_multiple_precision_check_correct_precision_and_value)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a("1.5",30), b("2.0");
	boost::multiprecision::mpfr_float expected_result("0.75");
	bertini::Float c = a/b;
	
	BOOST_CHECK_EQUAL(c.Precision(),30);
	BOOST_CHECK_EQUAL(*c.number_as_multiple_precision_,expected_result);
}








/****************
 **
 **  comparitors
 **
 ****************/




//      <     //
BOOST_AUTO_TEST_CASE(Float_less_than_double)
{
	bertini::Float a(1), b(1);
	
	BOOST_CHECK(!(a<b));
	BOOST_CHECK(!(b<a));
	
	BOOST_CHECK( (-a.MakeRandom() < b.MakeRandom()) );
}


BOOST_AUTO_TEST_CASE(Float_less_than_multiple_precision)
{
	bertini::Float a("1.0"), b("1.0");
	
	BOOST_CHECK(!(a<b));
	BOOST_CHECK(!(b<a));
	
	BOOST_CHECK( -a.MakeRandom() < b.MakeRandom() );
}



BOOST_AUTO_TEST_CASE(Float_less_than_mixed)
{
	bertini::Float a(1.0), b("1.0");
	
	BOOST_CHECK(!(a<b));
	BOOST_CHECK(!(b<a));
	
	BOOST_CHECK( -a.MakeRandom() < b.MakeRandom() );
}




//      >     //
BOOST_AUTO_TEST_CASE(Float_greater_than_double)
{
	bertini::Float a(1), b(1);
	
	BOOST_CHECK(!(a>b));
	BOOST_CHECK(!(b>a));
	
	BOOST_CHECK( a.MakeRandom() > -b.MakeRandom() );
}


BOOST_AUTO_TEST_CASE(Float_greater_than_multiple_precision)
{
	bertini::Float a("1.0"), b("1.0");
	
	BOOST_CHECK(!(a>b));
	BOOST_CHECK(!(b>a));
	
	BOOST_CHECK( a.MakeRandom() > -b.MakeRandom() );
}

BOOST_AUTO_TEST_CASE(Float_greater_than_mixed)
{
	
	bertini::Float a(1.0), b("1.0");
	
	BOOST_CHECK(!(a>b));
	BOOST_CHECK(!(b>a));
	
	BOOST_CHECK( a.MakeRandom() > -b.MakeRandom() );
}






//      <=     //
BOOST_AUTO_TEST_CASE(Float_less_than_or_equal_to_double)
{
	bertini::Float a(1), b(1);
	
	BOOST_CHECK(a<=b);
	BOOST_CHECK(b<=a);
	
	BOOST_CHECK( -a.MakeRandom() <= b.MakeRandom() );
}


BOOST_AUTO_TEST_CASE(Float_less_than_or_equal_to_multiple_precision)
{
	bertini::Float a("1.0"), b("1.0");
	
	BOOST_CHECK(a<=b);
	BOOST_CHECK(b<=a);
	
	BOOST_CHECK( -a.MakeRandom() <= b.MakeRandom() );
}



BOOST_AUTO_TEST_CASE(Float_less_than_or_equal_to_mixed)
{
	bertini::Float a(1.0), b("1.0");
	
	BOOST_CHECK((a<=b));
	BOOST_CHECK((b<=a));
	
	BOOST_CHECK( -a.MakeRandom() <= b.MakeRandom() );
}


//      >=     //
BOOST_AUTO_TEST_CASE(Float_greater_than_or_equal_to_double)
{
	bertini::Float a(1), b(1);
	
	BOOST_CHECK(a>=b);
	BOOST_CHECK(b>=a);
	
	BOOST_CHECK( a.MakeRandom() >= -b.MakeRandom() );
}


BOOST_AUTO_TEST_CASE(Float_greater_than_or_equal_to_multiple_precision)
{
	bertini::Float a("1.0"), b("1.0");
	
	BOOST_CHECK(a>=b);
	BOOST_CHECK(b>=a);
	
	BOOST_CHECK( a.MakeRandom() >= -b.MakeRandom() );
}



BOOST_AUTO_TEST_CASE(Float_greater_than_or_equal_to_mixed)
{
	bertini::Float a(1.0), b("1.0");
	
	BOOST_CHECK((a>=b));
	BOOST_CHECK((b>=a));
	
	BOOST_CHECK( a.MakeRandom() >= -b.MakeRandom() );
}










/***********************
 **
 **    the non-basic functions, such as square root, etc
 **
 *********************/










BOOST_AUTO_TEST_CASE(Float_square_root_double)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	bertini::Float a(2.0);
	
	BOOST_CHECK_EQUAL( sqrt(a),boost::multiprecision::mpfr_float("1.41421356237309504880168872421") );
}




BOOST_AUTO_TEST_CASE(Float_square_root_multiple_precision){
	
	
	boost::multiprecision::mpfr_float::default_precision(30);
	
	boost::multiprecision::mpfr_float expected_result( "1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641573");
	expected_result.precision(30);
	bertini::Float b("2.0",30);
	
	bertini::Float d = sqrt(b);
	BOOST_CHECK_EQUAL(d.Precision(),30);
	BOOST_CHECK_EQUAL(d.number_as_multiple_precision_->precision(),30);
	
	BOOST_CHECK_EQUAL(d, expected_result);
	
	
	boost::multiprecision::mpfr_float::default_precision(50);
	BOOST_CHECK_THROW(d = sqrt(b),std::exception);
}





BOOST_AUTO_TEST_CASE(Float_cosine_of_pi_is_minus_one)
{
	boost::multiprecision::mpfr_float::default_precision(30);
	
	BOOST_CHECK_EQUAL(cos(bertini::Float::pi<boost::multiprecision::mpfr_float>()) ,boost::multiprecision::mpfr_float(-1));
	
	boost::multiprecision::mpfr_float::default_precision(16);
	BOOST_CHECK_EQUAL(cos(bertini::Float::pi<double>()) ,double(-1));
}






BOOST_AUTO_TEST_CASE(make_pi_to_double_precision)
{
	
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float pi = bertini::Float::pi<double>();
	BOOST_CHECK_EQUAL( pi,3.141592653589793238462643383279502884197169399375105820);//89
}


BOOST_AUTO_TEST_CASE(make_pi_to_high_precision)
{
	boost::multiprecision::mpfr_float::default_precision(1000);
	bertini::Float pi = bertini::Float::pi<boost::multiprecision::mpfr_float>();
	
	BOOST_CHECK_EQUAL( pi,boost::multiprecision::mpfr_float("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420199"));//89
}



BOOST_AUTO_TEST_CASE(make_e_to_double_precision)
{
	boost::multiprecision::mpfr_float::default_precision(16);
	
	bertini::Float e = bertini::Float::e<double>();
	BOOST_CHECK_EQUAL(e,2.718281828459045235360287471352662497757247093699959574966967627724076630353);
}

BOOST_AUTO_TEST_CASE(make_e_to_high_precision)
{
	boost::multiprecision::mpfr_float::default_precision(1000);
	bertini::Float e = bertini::Float::e<boost::multiprecision::mpfr_float>();
	
	BOOST_CHECK_EQUAL(e,boost::multiprecision::mpfr_float("2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035159888851934580727386673858942287922849989208680582574927961048419844436346324496848756023362482704197862320900216099023530436994184914631409343173814364054625315209618369088870701676839642437814059271456354906130310720851038375051011574770417189861068739696552126715468895703503540"));
}





BOOST_AUTO_TEST_SUITE_END()










