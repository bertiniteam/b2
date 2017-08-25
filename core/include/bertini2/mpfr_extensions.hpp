//This file is part of Bertini 2.
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
//along with mpfr_extensions.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file mpfr_extensions.hpp 

\brief Extensions to the Boost.Multiprecision library.

Particularly includes Boost.Serialize code for the mpfr_float, gmp_rational, and gmp_int types.
*/

#ifndef BERTINI_MPFR_EXTENSIONS_HPP
#define BERTINI_MPFR_EXTENSIONS_HPP

#include "bertini2/config.h"

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#ifdef B2_FORBID_MIXED_ARITHMETIC
	#include "bertini2/forbid_double.hpp"
#endif

#include <boost/random.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>

#include "bertini2/double_extensions.hpp"


namespace bertini{

	

#ifdef BMP_EXPRESSION_TEMPLATES
	using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>; 
	using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_on>;
	using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_on>;
#else
	using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>; 
	using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_off>;
	using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_off>;
#endif

	inline auto DefaultPrecision()
	{
		return mpfr_float::default_precision();
	}

	inline void DefaultPrecision(unsigned prec)
	{
		mpfr_float::default_precision(prec);
	}

	/** 
	\brief Get the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline
	unsigned Precision(mpfr_float const& num)
	{
		return num.precision();
	}

	/** 
	\brief Change the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline void Precision(mpfr_float & num, unsigned prec)
	{
		num.precision(prec);
	}
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

#define EIGEN_DEVICE_FUNC // to make Eigen 3.3 happy... ugh, this is likely to break CUDA usage with Bertini2, if that ever happens.
#include <Eigen/src/Core/util/Macros.h>

#ifdef BMP_EXPRESSION_TEMPLATES
	#if (!EIGEN_VERSION_AT_LEAST(3,2,92)) // version of 3.3-beta1 is 3,2,92.
		namespace std{ 
using boost::multiprecision::min; //error receiver: please see https://svn.boost.org/trac/boost/ticket/11149 for information about these using statements in std namespace.  
using boost::multiprecision::max; //3 options: ./configure --disable-expression_templates, use Boost 1.61, or patch earlier Boost versions to resolve this.
		}
	#endif
#endif

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



#endif




