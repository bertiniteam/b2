//This file is part of Bertini 2.
//
//num_traits.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//num_traits.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with num_traits.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame

/**
\file num_traits.hpp 

\brief Provides an Eigen-like NumTraits struct for querying traits of a number type.  

The bertini::NumTraits struct provides NumDigits and NumFuzzyDigits functions.
*/

#ifndef BERTINI_NUM_TRAITS_HPP
#define BERTINI_NUM_TRAITS_HPP

#include <random>
#include <complex>
#include <cmath>
#include "bertini2/mpfr_complex.hpp"
#include "bertini2/mpfr_extensions.hpp"



namespace bertini
{

	using dbl = std::complex<double>;
	using mpfr = bertini::complex;

	template<typename T>
	T RandomUnit();

	template<typename T>
	struct NumTraits
	{};



	template <> struct NumTraits<double> 
	{
		inline static unsigned NumDigits()
		{
			return 16;
		}

		inline static unsigned NumFuzzyDigits()
		{
			return 14;
		}

		inline
		static unsigned TolToDigits(double tol)
		{
			return ceil(-log10(tol));
		}
	};


	template <> struct NumTraits<std::complex<double> > 
	{
		inline static unsigned NumDigits()
		{
			return 16;
		}

		inline static unsigned NumFuzzyDigits()
		{
			return 14;
		}
	};

	/**
	\brief Get the precision of a number.

	For doubles, this is trivially 16.
	*/
	inline
	unsigned Precision(const double num)
	{
		return 16;
	}

	/**
	\brief Get the precision of a number.

	For complex doubles, this is trivially 16.
	*/
	inline
	unsigned Precision(const std::complex<double> num)
	{
		return 16;
	}


	inline
	unsigned PrecisionIncrement()
	{
		return 10;
	}

	inline
	unsigned DoublePrecision()
	{
		return 16;
	}

	inline
	unsigned LowestMultiplePrecision()
	{
		return 20;
	}
		
	inline
	unsigned MaxPrecisionAllowed()
	{
		return 1000;
	}
	

	inline
	std::complex<double> rand_complex()
	{
		using std::abs;
		using std::sqrt;
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0,1.0);
		std::complex<double> returnme(distribution(generator), distribution(generator));
		return returnme / sqrt( abs(returnme));
	}

	template <> inline
	std::complex<double> RandomUnit<std::complex<double> >()
	{
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0,1.0);
		std::complex<double> returnme(distribution(generator), distribution(generator));
		return returnme / abs(returnme);
	}

	template <> 
	inline 
	bertini::complex RandomUnit<bertini::complex>()
	{
		return bertini::complex::RandomUnit();
	}
}// re: namespace bertini












namespace bertini {

	/** 
	\brief Get the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline
	unsigned Precision(mpfr_float const& num)
	{
		return num.precision();
	}
	
	
	
	template <> struct NumTraits<mpfr_float> 
	{
		inline static unsigned NumDigits()
		{
			return mpfr_float::default_precision();
		}

		inline static unsigned NumFuzzyDigits()
		{
			return mpfr_float::default_precision()-3;
		}

		inline
		static unsigned TolToDigits(mpfr_float tol)
		{
			auto b = ceil(-log10(tol));
			return b.convert_to<unsigned int>();
		}
	};	



	template <> struct NumTraits<bertini::complex> 
	{
		inline static unsigned NumDigits()
		{
			return mpfr_float::default_precision();
		}
	};
}

#endif


