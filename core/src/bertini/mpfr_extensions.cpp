//This file is part of Bertini 2.0.
//
//mpfr_extensions.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mpfr_extensions.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mpfr_extensions.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
// mpfr_extensions.cpp:  main source file for computational core of Bertini2


#include "mpfr_extensions.hpp"


namespace bertini {
	
	using mpfr_float = boost::multiprecision::mpfr_float;
	using mpq_rational = boost::multiprecision::mpq_rational;
	
	
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformUnitInterval()
	{
		static boost::uniform_01<T> uf;
		static boost::random::independent_bits_engine<
		boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
		> gen;
		return uf(gen);
	}
	
	template <unsigned int length_in_digits> mpfr_float RandomMpUniformUnitInterval();
	template <unsigned int length_in_digits> mpq_rational RandomMpUniformUnitInterval();
	
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformInInterval(const T & a, const T & b)
	{
		static boost::uniform_01<T> uf;
		static boost::random::independent_bits_engine<
		boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
		> gen;
		return (b-a)*uf(gen) + a;

	}
		
	template <unsigned int length_in_digits> mpfr_float RandomMpUniformInInterval();
	template <unsigned int length_in_digits> mpq_rational RandomMpUniformInInterval();


	template<typename T>
	T RandomMp(const T & a, const T & b)
	{
		auto num_digits = mpfr_float::default_precision() + 3;
		
		if (num_digits<=50)
			return RandomMpUniformInInterval<T,50>(a,b);
		else if (num_digits<=100)
			return RandomMpUniformInInterval<T,100>(a,b);
		else if (num_digits<=200)
			return RandomMpUniformInInterval<T,200>(a,b);
		else if (num_digits<=400)
			return RandomMpUniformInInterval<T,400>(a,b);
		else if (num_digits<=800)
			return RandomMpUniformInInterval<T,800>(a,b);
		else if (num_digits<=1600)
			return RandomMpUniformInInterval<T,1600>(a,b);
		else if (num_digits<=3200)
			return RandomMpUniformInInterval<T,3200>(a,b);
		else if (num_digits<=6400)
			return RandomMpUniformInInterval<T,6400>(a,b);
		else if (num_digits<=8000)
			return RandomMpUniformInInterval<T,8000>(a,b);
		else if (num_digits<=10000)
			return RandomMpUniformInInterval<T,10000>(a,b);
		else if (num_digits<=12000)
			return RandomMpUniformInInterval<T,12000>(a,b);
		else if (num_digits<=14000)
			return RandomMpUniformInInterval<T,14000>(a,b);
		else if (num_digits<=16000)
			return RandomMpUniformInInterval<T,16000>(a,b);
		else if (num_digits<=18000)
			return RandomMpUniformInInterval<T,18000>(a,b);
		else if (num_digits<=20000)
			return RandomMpUniformInInterval<T,20000>(a,b);
		else if (num_digits<=40000)
			return RandomMpUniformInInterval<T,40000>(a,b);
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000");
	}
		
	template mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b);
	template mpq_rational RandomMp(const mpq_rational & a, const mpq_rational & b);
	

	template <typename T>
	T RandomMp()
	{
		auto num_digits = mpfr_float::default_precision() + 3;
		
		if (num_digits<=50)
			return RandomMpUniformUnitInterval<T,50>();
		else if (num_digits<=100)
			return RandomMpUniformUnitInterval<T,100>();
		else if (num_digits<=200)
			return RandomMpUniformUnitInterval<T,200>();
		else if (num_digits<=400)
			return RandomMpUniformUnitInterval<T,400>();
		else if (num_digits<=800)
			return RandomMpUniformUnitInterval<T,800>();
		else if (num_digits<=1600)
			return RandomMpUniformUnitInterval<T,1600>();
		else if (num_digits<=3200)
			return RandomMpUniformUnitInterval<T,3200>();
		else if (num_digits<=6400)
			return RandomMpUniformUnitInterval<T,6400>();
		else if (num_digits<=8000)
			return RandomMpUniformUnitInterval<T,8000>();
		else if (num_digits<=10000)
			return RandomMpUniformUnitInterval<T,10000>();
		else if (num_digits<=12000)
			return RandomMpUniformUnitInterval<T,12000>();
		else if (num_digits<=14000)
			return RandomMpUniformUnitInterval<T,14000>();
		else if (num_digits<=16000)
			return RandomMpUniformUnitInterval<T,16000>();
		else if (num_digits<=18000)
			return RandomMpUniformUnitInterval<T,18000>();
		else if (num_digits<=20000)
			return RandomMpUniformUnitInterval<T,20000>();
		else if (num_digits<=40000)
			return RandomMpUniformUnitInterval<T,40000>();
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000");
	}

	template mpfr_float RandomMp();
	template mpq_rational RandomMp();
}




