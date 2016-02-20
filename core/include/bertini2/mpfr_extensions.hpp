//This file is part of Bertini 2.0.
//
//mpfr_extensions.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mpfr_extensions.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with .  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//  This file contains extensions to types provided by Boost.Multiprecision.
//
//

/**
\file mpfr_extensions.hpp 

\brief Extensions to the Boost.Multiprecision library.

Particularly includes Boost.Serialize code for the mpfr_float, gmp_rational, and gmp_int types.
*/

#ifndef BERTINI_MPFR_EXTENSIONS_HPP
#define BERTINI_MPFR_EXTENSIONS_HPP


#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/random.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>




namespace bertini{

	using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>;

	using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_off>;

	using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_off>;

	// using mpz_int = boost::multiprecision::mpz_int;
	// using mpq_rational = boost::multiprecision::mpq_rational;
}

// the following code block extends serialization to the mpfr_float class from boost::multiprecision
namespace boost { namespace serialization {
	/**
	 Save a mpfr_float type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0> const& r, unsigned /*version*/)
	{
		unsigned num_digits(r.precision());
		ar & num_digits;
		std::string tmp = r.str(0,std::ios::scientific);
		ar & tmp;
	}
	
	/**
	 Load a mpfr_float type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0>& r, unsigned /*version*/)
	{
		unsigned num_digits;
		ar & num_digits;
		r.precision(num_digits);
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}
	

	/**
	 Save a gmp_rational type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::gmp_rational const& r, unsigned /*version*/)
	{
		std::string tmp = r.str(0,std::ios::scientific);
		ar & tmp;
	}
	
	/**
	 Load a gmp_rational type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::gmp_rational& r, unsigned /*version*/)
	{
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}


	/**
	 Save a gmp_int type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::gmp_int const& r, unsigned /*version*/)
	{
		std::string tmp = r.str(0,std::ios::scientific);
		ar & tmp;
	}
	
	/**
	 Load a gmp_int type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::gmp_int& r, unsigned /*version*/)
	{
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}


} } // re: namespaces

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::mpfr_float_backend<0>)

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::gmp_rational)

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::gmp_int)


// if you wish to use et_on with Boost.Multiprecision with Eigen 3.2.x or earlier, you must apply a patch to Boost.MP related to bug 11149, and use the following non-standard lines.
// https://svn.boost.org/trac/boost/ticket/11149

// namespace std{ using boost::multiprecision::min; using
//  boost::multiprecision::max;
// }

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
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMp()
	{	

		using namespace boost::multiprecision;
   		using namespace boost::random;

   		static uniform_real_distribution<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<length_in_digits>, boost::multiprecision::et_off> > ur(0,1);
   		static independent_bits_engine<mt19937, length_in_digits*1000L/301L, mpz_int> gen;

		return mpfr_float(ur(gen));
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
	inline
	mpfr_float RandomMp()
	{
		auto num_digits = mpfr_float::default_precision() + 3;
	
		if (num_digits<=50)
			return RandomMp<50>();
		else if (num_digits<=100)
			return RandomMp<100>();
		else if (num_digits<=200)
			return RandomMp<200>();
		else if (num_digits<=400)
			return RandomMp<400>();
		else if (num_digits<=800)
			return RandomMp<800>();
		else if (num_digits<=1600)
			return RandomMp<1600>();
		else if (num_digits<=3200)
			return RandomMp<3200>();
		else if (num_digits<=6400)
			return RandomMp<6400>();
		else if (num_digits<=8000)
			return RandomMp<8000>();
		else if (num_digits<=10000)
			return RandomMp<10000>();
		else if (num_digits<=12000)
			return RandomMp<12000>();
		else if (num_digits<=14000)
			return RandomMp<14000>();
		else if (num_digits<=16000)
			return RandomMp<16000>();
		else if (num_digits<=18000)
			return RandomMp<18000>();
		else if (num_digits<=20000)
			return RandomMp<20000>();
		else if (num_digits<=40000)
			return RandomMp<40000>();
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.  If you have a better solution to this problem, please write the authors of this software.");
	}

	
	/**
	 \brief create a random number in a given interval, at the current default precision
	*/
	inline
	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b)
	{
		auto num_digits = mpfr_float::default_precision() + 3;
	 
		if (num_digits<=50)
			return RandomMp<50>(a,b);
		else if (num_digits<=100)
			return RandomMp<100>(a,b);
		else if (num_digits<=200)
			return RandomMp<200>(a,b);
		else if (num_digits<=400)
			return RandomMp<400>(a,b);
		else if (num_digits<=800)
			return RandomMp<800>(a,b);
		else if (num_digits<=1600)
			return RandomMp<1600>(a,b);
		else if (num_digits<=3200)
			return RandomMp<3200>(a,b);
		else if (num_digits<=6400)
			return RandomMp<6400>(a,b);
		else if (num_digits<=8000)
			return RandomMp<8000>(a,b);
		else if (num_digits<=10000)
			return RandomMp<10000>(a,b);
		else if (num_digits<=12000)
			return RandomMp<12000>(a,b);
		else if (num_digits<=14000)
			return RandomMp<14000>(a,b);
		else if (num_digits<=16000)
			return RandomMp<16000>(a,b);
		else if (num_digits<=18000)
			return RandomMp<18000>(a,b);
		else if (num_digits<=20000)
			return RandomMp<20000>(a,b);
		else if (num_digits<=40000)
			return RandomMp<40000>(a,b);
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.  If you have a better solution to this problem, please write the authors of this software.");
	}
	
	
	
	
	
} // re: namespace bertini



#endif




