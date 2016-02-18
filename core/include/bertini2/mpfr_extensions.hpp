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

#include <random>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>




namespace bertini{

	using mpfr_float = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_on>;

	using mpz_int = boost::multiprecision::mpz_int;
	using mpq_rational = boost::multiprecision::mpq_rational;
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




// from John Maddock, a boost multiprecision author, via email 2016.02.16 and 
// https://svn.boost.org/trac/boost/ticket/11149
namespace boost { namespace multiprecision {

template <class Backend, class tag, class A1, class A2, class A3, class A4> 
	inline number<Backend, et_on> min(const number<Backend, et_on>& arg, const detail::expression<tag, A1, A2, A3, A4>& a)
	{
		number<Backend, et_on> t(a);
		return (std::min)(arg, t);
	}
template <class tag, class A1, class A2, class A3, class A4, class Backend> 
	inline number<Backend, et_on> min(const detail::expression<tag, A1, A2, A3, A4>& arg, const number<Backend, et_on>& a)
	{
		number<Backend, et_on> t(arg);
		return (std::min)(arg, a);
	}
template <class tag, class A1, class A2, class A3, class A4, class tagb, class A1b, class A2b, class A3b, class A4b> 
	inline typename detail::expression<tag, A1, A2, A3, A4>::result_type min(const detail::expression<tag, A1, A2, A3, A4>& arg, const detail::expression<tagb, A1b, A2b, A3b, A4b>& a)
	{
		using N = typename detail::expression<tag, A1, A2, A3, A4>::result_type;
		N t1(arg), t2(a);
		return (std::min)(arg, a);
	}

template <class Backend, class tag, class A1, class A2, class A3, class A4>
	inline number<Backend, et_on> max(const number<Backend, et_on>& arg, const detail::expression<tag, A1, A2, A3, A4>& a)
	{
		number<Backend, et_on> t(a);
		return (std::max)(arg, t);
	}
template <class tag, class A1, class A2, class A3, class A4, class Backend>
	inline number<Backend, et_on> max(const detail::expression<tag, A1, A2, A3, A4>& arg, const number<Backend, et_on>& a)
	{
		number<Backend, et_on> t(arg);
		return (std::max)(arg, a);
	}
template <class tag, class A1, class A2, class A3, class A4, class tagb, class A1b, class A2b, class A3b, class A4b>
	inline typename detail::expression<tag, A1, A2, A3, A4>::result_type max(const detail::expression<tag, A1, A2, A3, A4>& arg, const detail::expression<tagb, A1b, A2b, A3b, A4b>& a)
	{	
		using N = typename detail::expression<tag, A1, A2, A3, A4>::result_type;
		N t1(arg), t2(a);
		return (std::max)(arg, a);
	}

} }

namespace bertini
{
		/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam T the number type to generate.
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 */
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformUnitInterval()
	{
		static boost::uniform_01<T> uf;
		static boost::random::independent_bits_engine<
			boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
														> gen;
		return uf(gen);
	}
	
	
	
	
	
	
	
	
	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam T the number type to generate.
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 \param a The left bound.
	 \param b The right bound.
	 */
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformInInterval(const T & a, const T & b)
	{
		static boost::uniform_01<T> uf;
		static boost::random::independent_bits_engine<
			boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int
														> gen;
		return (b-a)*uf(gen) + a;

	}
	
	
	
	/**
	 \brief create a random number in a given interval, at the current default precision
	 
	 \note this function calls the templated function RandomMpfrUniformInInterval.

	 \tparam T the number type to generate.
	 \param number_to_make_random The number whose contents you are overwriting with a random number.
	 */
	template <typename T>
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
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.");
	}
	
	
	/**
	 \brief create a random number, at the current default precision

	 \param number_to_make_random The number whose contents you are overwriting with a random number.
	 */
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
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.");
	}
	
	
} // re: namespace bertini



#endif




