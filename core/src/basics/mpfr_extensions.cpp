//This file is part of Bertini 2.
//
//mpfr_extensions.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mpfr_extensions.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mpfr_extensions.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file mpfr_extensions.cpp 

\brief Extensions to the Boost.Multiprecision library.

Particularly includes Boost.Serialize code for the mpfr_float, gmp_rational, and gmp_int types.
*/

#include "bertini2/mpfr_extensions.hpp"

namespace bertini {

	mpfr_float RandomMp()
	{
		auto num_digits = DefaultPrecision();

		if (num_digits<=50)
			return RandomMp<50>();
		else if (num_digits<=100)
			return RandomMp<100>();
		else if (num_digits<=200)
			return RandomMp<200>();
		else if (num_digits<=400)
			return RandomMp<400>();
		else if (num_digits<=800)
			return RandomMp<800>();
		else if (num_digits<=1600)
			return RandomMp<1600>();
		else if (num_digits<=3200)
			return RandomMp<3200>();
		else if (num_digits<=6400)
			return RandomMp<6400>();
		else if (num_digits<=8000)
			return RandomMp<8000>();
		else if (num_digits<=10000)
			return RandomMp<10000>();
		else if (num_digits<=12000)
			return RandomMp<12000>();
		else if (num_digits<=14000)
			return RandomMp<14000>();
		else if (num_digits<=16000)
			return RandomMp<16000>();
		else if (num_digits<=18000)
			return RandomMp<18000>();
		else if (num_digits<=20000)
			return RandomMp<20000>();
		else if (num_digits<=40000)
			return RandomMp<40000>();
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.  If you have a better solution to this problem, please write the authors of this software.");
	}




	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b)
	{
		auto num_digits = DefaultPrecision();

		if (num_digits<=50)
			return RandomMp<50>(a,b);
		else if (num_digits<=100)
			return RandomMp<100>(a,b);
		else if (num_digits<=200)
			return RandomMp<200>(a,b);
		else if (num_digits<=400)
			return RandomMp<400>(a,b);
		else if (num_digits<=800)
			return RandomMp<800>(a,b);
		else if (num_digits<=1600)
			return RandomMp<1600>(a,b);
		else if (num_digits<=3200)
			return RandomMp<3200>(a,b);
		else if (num_digits<=6400)
			return RandomMp<6400>(a,b);
		else if (num_digits<=8000)
			return RandomMp<8000>(a,b);
		else if (num_digits<=10000)
			return RandomMp<10000>(a,b);
		else if (num_digits<=12000)
			return RandomMp<12000>(a,b);
		else if (num_digits<=14000)
			return RandomMp<14000>(a,b);
		else if (num_digits<=16000)
			return RandomMp<16000>(a,b);
		else if (num_digits<=18000)
			return RandomMp<18000>(a,b);
		else if (num_digits<=20000)
			return RandomMp<20000>(a,b);
		else if (num_digits<=40000)
			return RandomMp<40000>(a,b);
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.  If you have a better solution to this problem, please write the authors of this software.");
	}







	void RandomMp(mpfr_float & a, unsigned num_digits)
	{
		if (num_digits<=50)
			RandomMp<50>(a);
		else if (num_digits<=100)
			RandomMp<100>(a);
		else if (num_digits<=200)
			RandomMp<200>(a);
		else if (num_digits<=400)
			RandomMp<400>(a);
		else if (num_digits<=800)
			RandomMp<800>(a);
		else if (num_digits<=1600)
			RandomMp<1600>(a);
		else if (num_digits<=3200)
			RandomMp<3200>(a);
		else if (num_digits<=6400)
			RandomMp<6400>(a);
		else if (num_digits<=8000)
			RandomMp<8000>(a);
		else if (num_digits<=10000)
			RandomMp<10000>(a);
		else if (num_digits<=12000)
			RandomMp<12000>(a);
		else if (num_digits<=14000)
			RandomMp<14000>(a);
		else if (num_digits<=16000)
			RandomMp<16000>(a);
		else if (num_digits<=18000)
			RandomMp<18000>(a);
		else if (num_digits<=20000)
			RandomMp<20000>(a);
		else if (num_digits<=40000)
			RandomMp<40000>(a);
		else
			throw std::out_of_range("requesting random long number of number of digits higher than 40000.  this can be remedied by adding more cases to the generating function RandomMp.  If you have a better solution to this problem, please write the authors of this software.");

		a.precision(num_digits);
	}


} // namespace bertini
