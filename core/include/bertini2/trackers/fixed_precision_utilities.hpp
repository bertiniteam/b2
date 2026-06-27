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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Tim Hodges, Colorado State University

#pragma once

#include "bertini2/trackers/config.hpp"


namespace bertini{ namespace tracking {
namespace fixed{

inline
void SetPrecision(SampCont<complex_mp> & samples, unsigned prec)
{
	for (auto& s : samples)
		for (unsigned ii=0; ii<s.size(); ii++)
			s(ii).precision(prec);
}

inline
void SetPrecision(TimeCont<complex_mp> & times, unsigned prec)
{
	for (auto& t : times)
		t.precision(prec);
}

inline
unsigned MaxPrecision(SampCont<complex_mp> const& samples)
{
	unsigned max_precision = 0;
	for (auto& s : samples)
		if(Precision(s(0)) > max_precision)
			max_precision = Precision(s(0));
	return max_precision;
}

inline
unsigned MaxPrecision(TimeCont<complex_mp> const& times)
{
	unsigned max_precision = 0;
	for (auto& t : times)
		if(Precision(t) > max_precision)
			max_precision = Precision(t);
	return max_precision;
}


//does not a thing, because cannot.
inline
unsigned EnsureAtUniformPrecision(TimeCont<complex_dbl> const & /*times*/, SampCont<complex_dbl> const & /*derivatives*/)
{
	return DoublePrecision();
}

inline
unsigned EnsureAtUniformPrecision(TimeCont<complex_dbl> const & /*times*/, SampCont<complex_dbl> const & /*samples*/, SampCont<complex_dbl> const & /*derivatives*/)
{
	return DoublePrecision(); 
}


//changes precision of complex_mp to highest needed precision for the samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<complex_mp> const & /*times*/, SampCont<complex_mp> const & samples)
{
	return MaxPrecision(samples);
}


//returns max precision of all samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<complex_mp> const & /*times*/, SampCont<complex_mp> const & samples, SampCont<complex_mp> const & derivatives)
{
	return max(MaxPrecision(samples),
	           MaxPrecision(derivatives)); 
}


}}} // re: namespaces


