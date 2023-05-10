//This file is part of Bertini 2.
//
//mpfr_complex.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mpfr_complex.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mpfr_complex.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file mpfr_complex.hpp

\brief The main multiprecision complex number type.  This is essentially boost::multiprecision' complex
*/


#ifndef BERTINI_MPFR_COMPLEX_HPP
#define BERTINI_MPFR_COMPLEX_HPP

#include "bertini2/config.h"



#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>
#include <assert.h>



#include "bertini2/mpfr_extensions.hpp"




#include <boost/multiprecision/mpc.hpp>

namespace bertini{

namespace bmp = boost::multiprecision;
using bmp::backends::mpc_complex_backend;

#ifdef BMP_EXPRESSION_TEMPLATES
	using mpfr_complex = bmp::number<mpc_complex_backend<0>, bmp::et_on >;
#else
	using mpfr_complex = bmp::number<mpc_complex_backend<0>, bmp::et_off >;
#endif

	inline auto DefaultPrecisionPolicy(){
		return bmp::variable_precision_options::preserve_source_precision;
	}

	

	// shamelessly adapted from the documentation for variable precision in Boost.Multiprecision.
	// see https://www.boost.org/doc/libs/1_82_0/libs/multiprecision/doc/html/boost_multiprecision/tut/variable.html
	struct scoped_mpfr_precision_options_this_thread
	{
	   boost::multiprecision::variable_precision_options saved_options;

	   scoped_mpfr_precision_options_this_thread(boost::multiprecision::variable_precision_options opts) : saved_options(mpfr_float::thread_default_variable_precision_options())
	   {
	      mpfr_float::thread_default_variable_precision_options(opts);
	   }

	   ~scoped_mpfr_precision_options_this_thread()
	   {
	      mpfr_float::thread_default_variable_precision_options(saved_options);
	   }

	   void reset(boost::multiprecision::variable_precision_options opts)
	   {
	      mpfr_float::thread_default_variable_precision_options(opts);
	   }

	};



	struct scoped_mpfr_precision_options_all_threads
	{
	   boost::multiprecision::variable_precision_options saved_options_all_threads;
	   boost::multiprecision::variable_precision_options saved_options_this_thread;

	   scoped_mpfr_precision_options_all_threads(boost::multiprecision::variable_precision_options opts) : 
	   		saved_options_all_threads(mpfr_float::default_variable_precision_options()), 
	   		saved_options_this_thread(mpfr_float::default_variable_precision_options())
	   {
	      mpfr_float::default_variable_precision_options(opts);
	      mpfr_float::thread_default_variable_precision_options(opts);
	   }

	   ~scoped_mpfr_precision_options_all_threads()
	   {
	      mpfr_float::default_variable_precision_options(saved_options_all_threads);
	      mpfr_float::thread_default_variable_precision_options(saved_options_this_thread);
	   }

	   void reset(boost::multiprecision::variable_precision_options opts)
	   {
	      mpfr_float::default_variable_precision_options(opts);
	      mpfr_float::thread_default_variable_precision_options(opts);
	   }

	};






	inline auto DefaultPrecision()
	{
		auto p = mpfr_float::default_precision();
		assert(p==mpfr_complex::default_precision() && "precision of real and complex multiprecision numbers have drifted...");
		return p;
	}

	inline void DefaultPrecision(unsigned prec)
	{
		mpfr_float::default_precision(prec);
		mpfr_complex::default_precision(prec);
	}


}


namespace boost { namespace serialization {
	/**
	 Save a mpc_complex type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::mpc_complex_backend<0> const& r, unsigned /*version*/)
	{
		unsigned num_digits(r.precision());
		ar & num_digits;
		std::string tmp = r.str(0,std::ios::scientific);
		ar & tmp;
	}

	/**
	 Load a mpc_complex type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::mpc_complex_backend<0>& r, unsigned /*version*/)
	{
		unsigned num_digits;
		ar & num_digits;
		r.precision(num_digits);
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}

}} // re: namespace boost::serialization


BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::mpc_complex_backend<0>);






namespace bertini{

	/** 
	\brief Get the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline
	auto Precision(mpfr_complex const& num)
	{
		return num.precision();
	}


	/** 
	\brief Change the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline void Precision(mpfr_complex & num, unsigned prec)
	{
		num.precision(prec);
	}

	inline
	bool isnan(mpfr_complex const& num){return isnan(num.real()) || isnan(num.imag());};



	using std::polar;
} // re: namespace bertini











#endif




