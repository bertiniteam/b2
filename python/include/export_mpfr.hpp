




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
		void ExportMpfr();
		
		using namespace boost::python;

		//typedef boost::multiprecision::number<boost::multiprecision::mpfr_float, boost::multiprecision::et_off> bmp;
//		typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0> , boost::multiprecision::et_off> bmp;
		typedef bertini::mpfr_float bmp;


		
	} //namespace python
} // namespace bertini

#endif

