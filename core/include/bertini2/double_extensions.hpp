//This file is part of Bertini 2.
//
//bertini2/double_extensions.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/double_extensions.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/double_extensions.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, University of Wisconsin Eau Claire

/**
\file bertini2/double_extensions.hpp 

\brief Provides bertini extensions to the double and complex<double> types.
*/

#ifndef BERTINI_DOUBLE_EXTENSIONS_HPP
#define BERTINI_DOUBLE_EXTENSIONS_HPP

#pragma once

#include <random>
#include <complex>


namespace bertini{

	/**
	\brief Overload * for unsigned * complex<double>
	*/
	inline
	std::complex<double> operator*(unsigned i, std::complex<double> z)
	{
		z*=i;
		return z;
	}

	/**
	an overload of isnan, for std::complex<double>
	*/
	inline
	bool isnan(std::complex<double> const& z)
	{
		using std::isnan;
		return isnan(z.real()) || isnan(z.imag());
	}


	/**
	\brief Gets you a random real number between -1 and 1, fwiw
	*/
	inline
	double RandReal()
	{
		static std::default_random_engine generator;
		static std::uniform_real_distribution<double> distribution(-1.0,1.0);
		return distribution(generator);
	}

	namespace{
		using dbl = std::complex<double>;
	}

	/**
	 Compute +,- integral powers of a std::complex<double> number.

	 This function recursively calls itself if the power is negative, by computing the power on the inverse.

	 \note This overload was removed from C++ in C++11, for some insane reason.  Here it is, back in black.
	 */
	inline dbl pow(const dbl & z, int power)
	{
		if (power < 0) {
			return pow(1./z, -power);
		}
		else if (power==0)
			return dbl(1,0);
		else if(power==1)
			return z;
		else if(power==2)
			return z*z;
		else if(power==3)
			return z*z*z;
		else
		{
			unsigned int p(power);
			dbl result(1,0), z_to_the_current_power_of_two = z;
			// have copy of p in memory, can freely modify it.
			do {
				if ( (p & 1) == 1 ) { // get the lowest bit of the number
					result *= z_to_the_current_power_of_two;
				}
				z_to_the_current_power_of_two *= z_to_the_current_power_of_two; // square z_to_the_current_power_of_two
			} while (p  >>= 1);
			
			return result;
		}
	}
	
} // namespace bertini


#endif // include guard

