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

#ifndef BERTINI_MPFR_EXTENSIONS_HPP
#define BERTINI_MPFR_EXTENSIONS_HPP


#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <random>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <eigen3/Eigen/Core>

#include <string>
#include <assert.h>




// the following code block extends serialization to the mpfr_float class from boost::multiprecision
namespace boost { namespace serialization {
	
	using mpfr_float = boost::multiprecision::mpfr_float;
	
	
	/**
	 Save a mpfr_float type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0> const& r, unsigned /*version*/)
	{
		std::string tmp = r.str(0,true);
		ar & tmp;
	}
	
	/**
	 Load a mpfr_float type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0>& r, unsigned /*version*/)
	{
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}
	
} } // re: namespaces

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::mpfr_float_backend<0>)






// reopen the Eigen namespace to inject this struct.
namespace Eigen {
	
	using boost::multiprecision::mpfr_float;
	using namespace boost::multiprecision;
	template<> struct NumTraits<boost::multiprecision::mpfr_float> : GenericNumTraits<boost::multiprecision::mpfr_float> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
		
		typedef boost::multiprecision::mpfr_float Real;
		typedef boost::multiprecision::mpfr_float NonInteger;
		typedef boost::multiprecision::mpfr_float Nested;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 20,
			AddCost = 30,
			MulCost = 40
		};
		
		
		inline static Real highest() {
			
			return (boost::multiprecision::mpfr_float(1) - epsilon()) * pow(boost::multiprecision::mpfr_float(2),mpfr_get_emax()-1);//);//boost::multiprecision::mpfr_float::default_precision());
		}
		
		inline static Real lowest() {
			return -highest();
		}
		
		inline static Real dummy_precision()
		{
			return pow( boost::multiprecision::mpfr_float(10),-int(boost::multiprecision::mpfr_float::default_precision()-3));
		}
		
		inline static Real epsilon()
		{
			return pow(boost::multiprecision::mpfr_float(10),-int(boost::multiprecision::mpfr_float::default_precision()));
		}
		//http://www.manpagez.com/info/mpfr/mpfr-2.3.2/mpfr_31.php
	};
}




namespace bertini {

	
	
	using mpfr_float = boost::multiprecision::mpfr_float;
	
	
	/**
	 \brief create a random number, at the current default precision
	 
	 \note this function calls the templated function RandomMpfr.
	 
	 \param number_to_make_random The number whose contents you are overwriting with a random number.
	 */
	template <typename T>
	T RandomMp();
	

	
	/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 */
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformUnitInterval();
	
	
	
	
	
	
	
	
	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 */
	template <typename T, unsigned int length_in_digits>
	T RandomMpUniformInInterval(const T & a, const T & b);
	
	
	
	/**
	 \brief create a random number in a given interval, at the current default precision
	 
	 \note this function calls the templated function RandomMpfrUniformInInterval.
	 
	 \param number_to_make_random The number whose contents you are overwriting with a random number.
	 */
	template <typename T>
	T RandomMp(const T & a, const T & b);
	
	
	
	
} // re: namespace bertini


#endif




