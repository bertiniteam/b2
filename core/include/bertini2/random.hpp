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
// Copyright(C) 2015 - 2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// danielle brake, university of wisconsin eau claire

/**
\file  bertini2/random.hpp 

\brief stuff for generating random numbers
*/

#ifndef BERTINI_RANDOM_HPP
#define BERTINI_RANDOM_HPP


#include <boost/random.hpp>

#include "bertini2/mpfr_complex.hpp"



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

   		static uniform_real_distribution<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<length_in_digits>, boost::multiprecision::et_off> > ur(0,1);
   		static independent_bits_engine<mt19937, length_in_digits*1000L/301L, mpz_int> gen;

		return ur(gen);
	}
	
	/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 */
	template <unsigned int length_in_digits>
	void RandomMp(mpfr_float & a)
	{	

		using namespace boost::multiprecision;
   		using namespace boost::random;

   		static uniform_real_distribution<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<length_in_digits>, boost::multiprecision::et_off> > ur(0,1);
   		static independent_bits_engine<mt19937, length_in_digits*1000L/301L, mpz_int> gen;

		a = ur(gen);
	}

	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 
	 \param a The left bound.
	 \param b The right bound.
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b)
	{
		return (b-a)*RandomMp<length_in_digits>()+a;
	}



	/**
	 \brief create a random number, at the current default precision
	 */
	mpfr_float RandomMp();

	/**
	 \brief create a random number in a given interval, at the current default precision
	*/
	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b);



	/**
	 \brief Set an existing mpfr_float to a random number, to a given precision.  

	 This function is how to get random numbers at a precision different from the current default.
	 */
	void RandomMp(mpfr_float & a, unsigned num_digits);

	

	
} // re: namespace bertini




namespace bertini{

namespace multiprecision{


using complex = bmp::number<mpc_complex_backend<0>, bmp::et_on >;
using bertini::RandomMp;




	inline 
	void RandomReal(complex & a, unsigned num_digits)
	{
		a.precision(num_digits);
		RandomMp(a.real(),num_digits);
		a.imag() = 0;
	}

	inline 
	void rand(complex & a, unsigned num_digits)
	{
		a.precision(num_digits);
		RandomMp(a.real(),num_digits);
		RandomMp(a.imag(),num_digits);
	}

	inline 
	void RandomComplex(complex & a, unsigned num_digits)
	{
		rand(a,num_digits);
	}


	/**
	 Produce a random complex number, to default precision.
	 */
	inline complex rand()
	{
		complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		returnme /= sqrt( abs(returnme));
		return returnme;
	}
	
	inline complex RandomUnit()
	{
		complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		returnme /= abs(returnme);
		return returnme;
	}
	/**
	 Produce a random real number \f$\in [-1,\,1]\f$, to current default precision. 
	 */
	inline complex RandomReal()
	{
		using std::sqrt;
		complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		returnme /= sqrt( abs(returnme));
		return returnme;
	}
		


	inline 
	complex RandomComplex(unsigned num_digits)
	{
		complex z;
		RandomComplex(z, num_digits);
		return z;
	}

	inline 
	void RandomUnit(complex & a, unsigned num_digits)
	{
		auto prev_precision = DefaultPrecision();

		a.precision(num_digits);
		RandomMp(a.real(),num_digits);
		RandomMp(a.imag(),num_digits);
		a /= abs(a);

		DefaultPrecision(prev_precision);
	}

	inline 
	complex RandomUnit(unsigned num_digits)
	{
		complex a;
		RandomUnit(a,num_digits);
		return a;
	}

}  // namespace multiprecision


} // namespaces
// }



#endif




