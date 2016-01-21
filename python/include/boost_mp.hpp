




#ifndef BERTINI_PYTHON_BOOST_MP_HPP
#define BERTINI_PYTHON_BOOST_MP_HPP



#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
// #include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <boost/python/wrapper.hpp>

#include <boost/python/operators.hpp>
#include <boost/operators.hpp>


#include <bertini2/mpfr_extensions.hpp>

#include <sstream>



namespace bertini{
	namespace python{
		void ExportMpfrFloat();
		
		using namespace boost::python;

		//typedef boost::multiprecision::number<boost::multiprecision::mpfr_float, boost::multiprecision::et_off> bmp;
		typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0> , boost::multiprecision::et_off> bmp;


		namespace detail{
			// a short little function wrapping the metaprogramming stuff from boost::multiprecision, so we can expose sqrt.  if you, reading this right now, know how to make these work with expression templates, please help the project by implementing it, and issuing a pull request.
			inline
			bmp bmpsqrt(bmp const& n)
			{
				return sqrt(n);
			}


			// inline
			// std::string bmp_format(bmp const& n, char* fmt)
			// {
			// 	std::stringstream ss;

			// 	unsigned digits = 5;




			// }


			// a free function getting a number as a string
			inline
			std::string bmp_as_str(bmp const& n)
			{
				std::stringstream ss;
				ss << n;
				return ss.str();
			}

			// a free function getting a number as a string, to full precision
			inline
			std::string full_precision_string(bmp const& n)
			{	
				std::stringstream ss;
				ss << n.str(0,std::ios::scientific);
				return ss.str();
			}
		}
		
	} //namespace python
} // namespace bertini

#endif

