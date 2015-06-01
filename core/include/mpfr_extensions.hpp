
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
	mpfr_float RandomMpfr();
	

	
	/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMpfrUniformUnitInterval();
	
	
	
	
	
	
	
	
	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \return number_to_make_random The number which you desire to populate with a random number.
	 
	 */
	template <unsigned int length_in_digits>
	mpfr_float RandomMpfrUniformInInterval(const mpfr_float & a, const mpfr_float & b);
	
	
	
	/**
	 \brief create a random number in a given interval, at the current default precision
	 
	 \note this function calls the templated function RandomMpfrUniformInInterval.
	 
	 \param number_to_make_random The number whose contents you are overwriting with a random number.
	 */
	mpfr_float RandomMpfr(const mpfr_float & a, const mpfr_float & b);
	
	
	
	
} // re: namespace bertini


#endif




