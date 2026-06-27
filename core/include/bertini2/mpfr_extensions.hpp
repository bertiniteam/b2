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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin - eau claire

/**
\file mpfr_extensions.hpp 

\brief Extensions to the Boost.Multiprecision library.

Particularly includes Boost.Serialize code for the real_mp, gmp_rational, and gmp_int types.
*/

#ifndef BERTINI_MPFR_EXTENSIONS_HPP
#define BERTINI_MPFR_EXTENSIONS_HPP



#include <boost/multiprecision/mpfr.hpp>

#ifdef B2_FORBID_MIXED_ARITHMETIC
	#include "bertini2/forbid_mixed_arithmetic.hpp"
#endif

#include <boost/random.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>

#include "bertini2/double_extensions.hpp"

namespace bertini{
#ifdef BMP_EXPRESSION_TEMPLATES
	using real_mp = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>; 
	using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_on>;
	using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_on>;
#else
	using real_mp = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>; 
	using mpz_int = boost::multiprecision::number<boost::multiprecision::backends::gmp_int, boost::multiprecision::et_off>;
	using mpq_rational = boost::multiprecision::number<boost::multiprecision::backends::gmp_rational, boost::multiprecision::et_off>;
#endif

	/** 
	\brief Get the precision of a real number.

	For mpfr_floats, this calls the precision member method for real_mp.
	*/
	inline
	auto Precision(real_mp const& num)
	{
		return num.precision();
	}

	/** 
	\brief Change the precision of a real number.

	For mpfr_floats, this calls the precision member method for real_mp.
	*/
	inline void Precision(real_mp & num, unsigned prec)
	{
		num.precision(prec);
	}

}


// overloads of the `max` function.
namespace bertini{
	/**
	\brief Three-argument form of `max`.

	\param a Input one
	\param b Input two
	\param c Input three
	*/
	template <typename T>
	T max(T const& a, T const& b, T const& c)
	{
		using std::max;
		return max(max(a,b),c);
	}

	/**
	\brief Four-argument form of `max`.

	\param a Input one
	\param b Input two
	\param c Input three
	\param d Input four
	*/
	template <typename T>
	T max(T const& a, T const& b, T const& c, T const& d)
	{
		using std::max;
		using bertini::max;
		return max(max(a,b,c),d);
	}

	/**
	\brief Five-argument form of `max`.

	\param a Input one
	\param b Input two
	\param c Input three
	\param d Input four
	\param e Input five
	*/
	template <typename T>
	T max(T const& a, T const& b, T const& c, T const& d, T const& e)
	{
		using std::max;
		using bertini::max;
		return max(max(a,b,c,d),e);
	}
}



// the following code block extends serialization to the real_mp class from boost::multiprecision
namespace boost { namespace serialization {
	/**
	 Save a real_mp type to a boost archive.
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
	 Load a real_mp type from a boost archive.
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


// Eigen 3.2.x workaround for Boost.Multiprecision expression templates
// removed -- project now requires Eigen >= 3.3.



#endif




