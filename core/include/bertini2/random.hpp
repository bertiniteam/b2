//This file is part of Bertini 2.
//
// bertini2/random.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// bertini2/random.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  bertini2/random.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file  bertini2/random.hpp 

\brief stuff for generating random numbers
*/

#ifndef BERTINI_RANDOM_HPP
#define BERTINI_RANDOM_HPP




#include "bertini2/mpfr_complex.hpp"
#include <boost/random.hpp>


namespace bertini
{
	/**
	Generate a random integer number between -10^digits and 10^digits
	*/
	template <unsigned long digits = 50>
	inline
	mpz_int RandomInt()
	{
		using namespace boost::random;
   		static mt19937 mt;
	    static uniform_int_distribution<mpz_int> ui(-(mpz_int(1) << digits*1000L/301L), mpz_int(1) << digits*1000L/301L);
	    return ui(mt);
	}
	
	
	/**
	Generate a random rational number with numerator and denomenator between -10^digits and 10^digits
	*/
	template <unsigned long digits = 50>
	mpq_rational RandomRat()
	{
   		using namespace boost::random;
   		static mt19937 mt;
	    static uniform_int_distribution<mpz_int> ui(-(mpz_int(1) << digits*1000L/301L), mpz_int(1) << digits*1000L/301L);
	    return mpq_rational(ui(mt),ui(mt));
	}


	/**
	 Produce a random number with at length_in_digits non-zero digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMp()
	{	

		using namespace boost::multiprecision;
   		using namespace boost::random;

   		static uniform_real_distribution<number<mpfr_float_backend<length_in_digits>, et_on> > distribution(0,1);
   		static independent_bits_engine<mt19937, length_in_digits*1000L/301L, mpz_int> bit_generator;

		mpfr_float a{distribution(bit_generator)};
		return a;
	}
	
	/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \param a the number which will be assigned in this call
	 */
	template <unsigned int length_in_digits>
	void RandomMpAssign(mpfr_float & a)
	{	
		a = RandomMp<length_in_digits>();
	}

	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 
	 \param a The left bound.
	 \param b The right bound.
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMp(const mpfr_float & left, const mpfr_float & right)
	{
		return (right-left)*RandomMp<length_in_digits>()+left;
	}



	/**
	 \brief create a random number, at the current default precision
	 */
	mpfr_float RandomMp();

	/**
	 \brief create a random number, at the specified precision

	 \param num_digits the precision that you desire.  
	 */
	mpfr_float RandomMp(unsigned num_digits);

	/**
	 \brief create a random number in a given interval, at the current default precision
	*/
	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b);

	/**
	 \brief create a random number in a given interval, at the specified precision
	*/
	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b, unsigned num_digits);

	/**
	 \brief Set an existing mpfr_float to a random number, to a given precision.  

	 This function is how to get random numbers at a precision different from the current default.
	 */
	void RandomMpAssign(mpfr_float & a, unsigned num_digits);

	

	
} // re: namespace bertini




namespace bertini{

namespace multiprecision{


using complex = bertini::mpfr_complex;
using bertini::RandomMp;



	/**
	 Assign to a random real number \f$\in [-1,\,1]\f$, to current default precision. 
	 */
	inline 
	void RandomRealAssign(complex & a, unsigned num_digits)
	{
		auto cached = DefaultPrecision();
		DefaultPrecision(num_digits);
		complex temp(RandomMp(mpfr_float(-1),mpfr_float(1),num_digits)); // ,0
		a.swap(temp);
		DefaultPrecision(cached);
	}

	/**
	 Produce a random real number \f$\in [-1,\,1]\f$, to current default precision. 
	 */
	inline complex RandomReal()
	{
		return complex(RandomMp(mpfr_float(-1),mpfr_float(1))); // ,0
	}
	
	/**
	 Produce a random real number \f$\in [-1,\,1]\f$, to specified precision. 
	 */
	inline complex RandomReal(unsigned num_digits)
	{
		auto cached = DefaultPrecision();
		DefaultPrecision(num_digits);
		auto result = complex(RandomMp(mpfr_float(-1),mpfr_float(1),num_digits));// ,0
		DefaultPrecision(cached);
		return result;
	}



	
		


	/**
	 Produce a random complex number, to default precision.
	 */
	inline complex rand()
	{
		return complex( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
	}


	/**
	 Produce a random unit complex number, to default precision.
	 */
	inline complex rand_unit()
	{
		complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		return returnme / sqrt( abs(returnme));
	}

	inline complex RandomUnit()
	{
		return rand_unit();
	}

	inline 
	void rand_assign(complex & a, unsigned num_digits)
	{
		auto cached = DefaultPrecision();
		DefaultPrecision(num_digits);
		
		mpfr_complex temp( RandomMp(num_digits), RandomMp(num_digits) );
		a = std::move(temp);
		DefaultPrecision(cached);
	}

	inline 
	void RandomComplexAssign(complex & a, unsigned num_digits)
	{
		rand_assign(a,num_digits);
	}

	inline 
	complex RandomComplex(unsigned num_digits)
	{
		complex z;
		RandomComplexAssign(z, num_digits);
		return z;
	}

	inline 
	void RandomUnitAssign(complex & a, unsigned num_digits)
	{
		auto cached = DefaultPrecision();
		DefaultPrecision(num_digits);
		a.precision(num_digits);
		
		complex temp(RandomMp(num_digits),RandomMp(num_digits));
		a = std::move(temp/sqrt(abs(temp)));
		DefaultPrecision(cached);
	}

	inline 
	complex RandomUnit(unsigned num_digits)
	{
		complex a;
		RandomUnitAssign(a,num_digits);
		return a;
	}

}  // namespace multiprecision


} // namespaces
// }



#endif




