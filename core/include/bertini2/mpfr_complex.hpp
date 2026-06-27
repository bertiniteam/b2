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
// Copyright(C) Bertini2 Development Team
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
	using complex_mp = bmp::number<mpc_complex_backend<0>, bmp::et_on >;
#else
	using complex_mp = bmp::number<mpc_complex_backend<0>, bmp::et_off >;
#endif

	inline auto DefaultPrecisionPolicy(){
		return bmp::variable_precision_options::preserve_related_precision;
	}




	// shamelessly adapted from the documentation for variable precision in Boost.Multiprecision.
	// see https://www.boost.org/doc/libs/1_82_0/libs/multiprecision/doc/html/boost_multiprecision/tut/variable.html
	struct scoped_mpfr_precision_options_this_thread
	{
	   boost::multiprecision::variable_precision_options saved_options;

	   scoped_mpfr_precision_options_this_thread(boost::multiprecision::variable_precision_options opts) : saved_options(real_mp::thread_default_variable_precision_options())
	   {
	      real_mp::thread_default_variable_precision_options(opts);
	   }

	   ~scoped_mpfr_precision_options_this_thread()
	   {
	      real_mp::thread_default_variable_precision_options(saved_options);
	   }

	   void reset(boost::multiprecision::variable_precision_options opts)
	   {
	      real_mp::thread_default_variable_precision_options(opts);
	   }

	};



	struct scoped_mpfr_precision_options_all_threads
	{
	   boost::multiprecision::variable_precision_options saved_options_all_threads;
	   boost::multiprecision::variable_precision_options saved_options_this_thread;

	   scoped_mpfr_precision_options_all_threads(boost::multiprecision::variable_precision_options opts) :
	   		saved_options_all_threads(real_mp::default_variable_precision_options()),
	   		saved_options_this_thread(real_mp::default_variable_precision_options())
	   {
	      real_mp::default_variable_precision_options(opts);
	      real_mp::thread_default_variable_precision_options(opts);
	   }

	   ~scoped_mpfr_precision_options_all_threads()
	   {
	      real_mp::default_variable_precision_options(saved_options_all_threads);
	      real_mp::thread_default_variable_precision_options(saved_options_this_thread);
	   }

	   void reset(boost::multiprecision::variable_precision_options opts)
	   {
	      real_mp::default_variable_precision_options(opts);
	      real_mp::thread_default_variable_precision_options(opts);
	   }

	};






	inline auto DefaultPrecision()
	{
		auto p = real_mp::default_precision();
		assert(p==complex_mp::default_precision() && "precision of real and complex multiprecision numbers have drifted...");
		return p;
	}

	inline void DefaultPrecision(unsigned prec)
	{
		real_mp::default_precision(prec);
		complex_mp::default_precision(prec);
		// Boost.Multiprecision >= 1.87 reads thread_default_precision when
		// default-constructing mpfr temporaries (including those created by
		// boost::python const& extractors). If left at 0, mpfr_init2 aborts.
		// Set both so that calls to default_precision() align thread-local
		// with the static default. See commit 3111255b for original context.
		real_mp::thread_default_precision(prec);
		complex_mp::thread_default_precision(prec);
#ifdef BMP_EXPRESSION_TEMPLATES
		// With ET on, the default global policy for mpc_complex_backend is preserve_related_precision
		// and for mpfr_float_backend it is preserve_target_precision. Both differ from what we want.
		// Setting preserve_related_precision per-thread mirrors the mpc_complex_backend global default
		// and ensures copy/construction semantics match the et_off behavior. Specifically:
		//   - mpc_complex copy ctor uses preserve_related_precision() (>= 3) to preserve source precision
		//   - assign_components_set_precision uses preserve_component_precision() (>= 2) to resize
		//     a complex from real components at higher-than-default precision
		//   - mpc_complex = real_mp uses preserve_component_precision() (>= 2) to resize
		// preserve_related_precision satisfies all three thresholds.
		// For real_mp, preserve_related_precision also enables source-precision-preserving copies.
		real_mp::thread_default_variable_precision_options(
			bmp::variable_precision_options::preserve_related_precision);
		complex_mp::thread_default_variable_precision_options(
			bmp::variable_precision_options::preserve_related_precision);
#endif
	}

	// Sets thread-local precision only — does NOT write the global default_precision.
	// Safe to call concurrently from multiple std::thread workers tracking at different
	// precisions. Use instead of DefaultPrecision() inside per-thread tracking loops.
	inline void SetThreadPrecision(unsigned prec)
	{
		real_mp::thread_default_precision(prec);
		complex_mp::thread_default_precision(prec);
	}

	inline unsigned ThreadPrecision()
	{
		return static_cast<unsigned>(real_mp::thread_default_precision());
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

	For mpfr_floats, this calls the precision member method for real_mp.
	*/
	inline
	auto Precision(complex_mp const& num)
	{
		return num.precision();
	}


	/**
	\brief Change the precision of a number.

	For mpfr_floats, this calls the precision member method for real_mp.
	*/
	inline void Precision(complex_mp & num, unsigned prec)
	{
		num.precision(prec);
	}

	inline
	bool isnan(complex_mp const& num){return isnan(num.real()) || isnan(num.imag());};



	using std::polar;
} // re: namespace bertini











#endif




