//This file is part of Bertini 2.
//
//random.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//random.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with random.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file random.cpp 

\brief stuff to make random numbers in bertini2
*/

#include "bertini2/random.hpp"

namespace bertini {

	mpfr_float RandomMp()
	{
		return RandomMp(bertini::DefaultPrecision());
	}

	mpfr_float RandomMp(unsigned num_digits)
	{
		
		mpfr_float a;
		if (num_digits<=50)
			a = std::move(RandomMp<50>());
		else if (num_digits<=100)
			a = std::move(RandomMp<100>());
		else if (num_digits<=200)
			a = std::move(RandomMp<200>());
		else if (num_digits<=400)
			a = std::move(RandomMp<400>());
		else if (num_digits<=800)
			a = std::move(RandomMp<800>());
		else if (num_digits<=1600)
			a = std::move(RandomMp<1600>());
		else if (num_digits<=3200)
			a = std::move(RandomMp<3200>());
		else if (num_digits<=6400)
			a = std::move(RandomMp<6400>());
		else if (num_digits<=8000)
			a = std::move(RandomMp<8000>());
		else if (num_digits<=10000)
			a = std::move(RandomMp<10000>());
		else if (num_digits<=12000)
			a = std::move(RandomMp<12000>());
		else if (num_digits<=14000)
			a = std::move(RandomMp<14000>());
		else if (num_digits<=16000)
			a = std::move(RandomMp<16000>());
		else if (num_digits<=18000)
			a = std::move(RandomMp<18000>());
		else if (num_digits<=20000)
			a = std::move(RandomMp<20000>());
		else if (num_digits<=40000)
			a = std::move(RandomMp<40000>());
		else
			throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
		a.precision(num_digits);
		return a;
	}


	void RandomMpAssign(mpfr_float & a, unsigned num_digits)
	{
		

		mpfr_float temp;
		temp = RandomMp(num_digits);
		a = std::move(temp);
	}




	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b)
	{
		return RandomMp(a,b,bertini::DefaultPrecision());
	}

	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b, unsigned num_digits)
	{
		
		mpfr_float result;
		if (num_digits<=50)
			result = std::move(RandomMp<50>(a,b));
		else if (num_digits<=100)
			result = std::move(RandomMp<100>(a,b));
		else if (num_digits<=200)
			result = std::move(RandomMp<200>(a,b));
		else if (num_digits<=400)
			result = std::move(RandomMp<400>(a,b));
		else if (num_digits<=800)
			result = std::move(RandomMp<800>(a,b));
		else if (num_digits<=1600)
			result = std::move(RandomMp<1600>(a,b));
		else if (num_digits<=3200)
			result = std::move(RandomMp<3200>(a,b));
		else if (num_digits<=6400)
			result = std::move(RandomMp<6400>(a,b));
		else if (num_digits<=8000)
			result = std::move(RandomMp<8000>(a,b));
		else if (num_digits<=10000)
			result = std::move(RandomMp<10000>(a,b));
		else if (num_digits<=12000)
			result = std::move(RandomMp<12000>(a,b));
		else if (num_digits<=14000)
			result = std::move(RandomMp<14000>(a,b));
		else if (num_digits<=16000)
			result = std::move(RandomMp<16000>(a,b));
		else if (num_digits<=18000)
			result = std::move(RandomMp<18000>(a,b));
		else if (num_digits<=20000)
			result = std::move(RandomMp<20000>(a,b));
		else if (num_digits<=40000)
			result = std::move(RandomMp<40000>(a,b));
		else
			throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
		result.precision(num_digits);
		return result;
	}






	


} // namespace bertini
