//This file is part of Bertini 2.
//
//fixed_precision_utilities.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_precision_utilities.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_precision_utilities.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University

#pragma once

#include "bertini2/trackers/config.hpp"


namespace bertini{ namespace tracking {
namespace fixed{

inline
void SetPrecision(SampCont<mpfr> & samples, unsigned prec)
{
	for (auto& s : samples)
		for (unsigned ii=0; ii<s.size(); ii++)
			s(ii).precision(prec);
}

inline
void SetPrecision(TimeCont<mpfr> & times, unsigned prec)
{
	for (auto& t : times)
		t.precision(prec);
}

inline
unsigned MaxPrecision(SampCont<mpfr> const& samples)
{
	unsigned max_precision = 0;
	for (auto& s : samples)
		if(Precision(s(0)) > max_precision)
			max_precision = Precision(s(0));
	return max_precision;
}

inline
unsigned MaxPrecision(TimeCont<mpfr> const& times)
{
	unsigned max_precision = 0;
	for (auto& t : times)
		if(Precision(t) > max_precision)
			max_precision = Precision(t);
	return max_precision;
}


//does not a thing, because cannot.
inline
unsigned EnsureAtUniformPrecision(TimeCont<dbl> const & times, SampCont<dbl> const & derivatives)
{
	return DoublePrecision();
}

inline
unsigned EnsureAtUniformPrecision(TimeCont<dbl> const & times, SampCont<dbl> const & samples, SampCont<dbl> const & derivatives)
{
	return DoublePrecision(); 
}


//changes precision of mpfr to highest needed precision for the samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> const & times, SampCont<mpfr> const & samples)
{
	return MaxPrecision(samples);
}


//returns max precision of all samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> const & times, SampCont<mpfr> const & samples, SampCont<mpfr> const & derivatives)
{
	return max(MaxPrecision(samples),
	           MaxPrecision(derivatives)); 
}


}}} // re: namespaces


