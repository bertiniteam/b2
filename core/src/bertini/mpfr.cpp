#include "mpfr_extensions.hpp"


namespace bertini {
	
	using mpfr_float = boost::multiprecision::mpfr_float;

	
	
	template <unsigned int length_in_digits>
	mpfr_float RandomMpfrUniformUnitInterval()
	{
		static boost::uniform_01<mpfr_float> uf;
		static boost::random::independent_bits_engine<
		boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
		> gen;
		return uf(gen);
	}
	
	
	
	template <unsigned int length_in_digits>
	mpfr_float RandomMpfrUniformInInterval(const mpfr_float & a, const mpfr_float & b)
	{
		static boost::uniform_01<mpfr_float> uf;
		static boost::random::independent_bits_engine<
		boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
		> gen;
		return (b-a)*uf(gen) + a;

	}
	
	
	mpfr_float RandomMpfr(const mpfr_float & a, const mpfr_float & b)
	{
		auto num_digits = mpfr_float::default_precision() + 3;
		
		if (num_digits<=50)
			return RandomMpfrUniformInInterval<50>(a,b);
		else if (num_digits<=100)
			return RandomMpfrUniformInInterval<100>(a,b);
		else if (num_digits<=200)
			return RandomMpfrUniformInInterval<200>(a,b);
		else if (num_digits<=400)
			return RandomMpfrUniformInInterval<400>(a,b);
		else if (num_digits<=800)
			return RandomMpfrUniformInInterval<800>(a,b);
		else if (num_digits<=1600)
			return RandomMpfrUniformInInterval<1600>(a,b);
		else if (num_digits<=3200)
			return RandomMpfrUniformInInterval<3200>(a,b);
		else if (num_digits<=6400)
			return RandomMpfrUniformInInterval<6400>(a,b);
		else if (num_digits<=8000)
			return RandomMpfrUniformInInterval<8000>(a,b);
		else if (num_digits<=10000)
			return RandomMpfrUniformInInterval<10000>(a,b);
		else if (num_digits<=12000)
			return RandomMpfrUniformInInterval<12000>(a,b);
		else if (num_digits<=14000)
			return RandomMpfrUniformInInterval<14000>(a,b);
		else if (num_digits<=16000)
			return RandomMpfrUniformInInterval<16000>(a,b);
		else if (num_digits<=18000)
			return RandomMpfrUniformInInterval<18000>(a,b);
		else if (num_digits<=20000)
			return RandomMpfrUniformInInterval<20000>(a,b);
		else if (num_digits<=40000)
			return RandomMpfrUniformInInterval<40000>(a,b);
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000");
	}
	
	
	
	mpfr_float RandomMpfr()
	{
		auto num_digits = mpfr_float::default_precision() + 3;
		
		if (num_digits<=50)
			return RandomMpfrUniformUnitInterval<50>();
		else if (num_digits<=100)
			return RandomMpfrUniformUnitInterval<100>();
		else if (num_digits<=200)
			return RandomMpfrUniformUnitInterval<200>();
		else if (num_digits<=400)
			return RandomMpfrUniformUnitInterval<400>();
		else if (num_digits<=800)
			return RandomMpfrUniformUnitInterval<800>();
		else if (num_digits<=1600)
			return RandomMpfrUniformUnitInterval<1600>();
		else if (num_digits<=3200)
			return RandomMpfrUniformUnitInterval<3200>();
		else if (num_digits<=6400)
			return RandomMpfrUniformUnitInterval<6400>();
		else if (num_digits<=8000)
			return RandomMpfrUniformUnitInterval<8000>();
		else if (num_digits<=10000)
			return RandomMpfrUniformUnitInterval<10000>();
		else if (num_digits<=12000)
			return RandomMpfrUniformUnitInterval<12000>();
		else if (num_digits<=14000)
			return RandomMpfrUniformUnitInterval<14000>();
		else if (num_digits<=16000)
			return RandomMpfrUniformUnitInterval<16000>();
		else if (num_digits<=18000)
			return RandomMpfrUniformUnitInterval<18000>();
		else if (num_digits<=20000)
			return RandomMpfrUniformUnitInterval<20000>();
		else if (num_digits<=40000)
			return RandomMpfrUniformUnitInterval<40000>();
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000");
	}

}




