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

/**
\file tracking/adaptive_precision_utilities.hpp 

\brief Functions for dealing with precisions of objects, particularly in the context of adaptive precision.

These definitions are provided for use in tracking and endgames.  This allows a uniform interface regardless of tracker details.
*/

#pragma once

#include "bertini2/trackers/config.hpp"


namespace bertini{ namespace tracking { namespace adaptive{

/**
\brief Sets the precision of each space sample to be of input precision.

\param samples The samples of which to change precision.
\param prec The new precision the samples should have.
*/
inline
void SetPrecision(SampCont<mpfr> & samples, unsigned prec)
{
	for (auto& s : samples)
		for (unsigned ii=0; ii<s.size(); ii++)
			s(ii).precision(prec);
}

/**
\brief Sets the precision of each time sample to be of input precision.

\param times The times of which to change precision.
\param prec The new precision the times should have.
*/
inline
void SetPrecision(TimeCont<mpfr> & times, unsigned prec)
{
	for (auto& t : times)
		t.precision(prec);
}

/**
\brief Get the maximum precision among an input set of space samples.  

This is computed based on the first coordinate.  Length-zero samples will cause either an assert or a throw (from Eigen, not Bertini2).

\param samples Some complex samples in space, obtained from tracking probably.
\return The maximum precision among those samples.
*/
inline
unsigned MaxPrecision(SampCont<mpfr> const& samples)
{
	unsigned max_precision = 0;
	for (const auto& s : samples)
		if(Precision(s(0)) > max_precision)
			max_precision = Precision(s(0));
	return max_precision;
}

/**
\brief Gets the maximum precision among a set of times.

\param times Some complex times.
\return The maximum precision among those times.
*/
inline
unsigned MaxPrecision(TimeCont<mpfr> const& times)
{
	unsigned max_precision = 0;
	for (const auto& t : times)
		if(Precision(t) > max_precision)
			max_precision = Precision(t);
	return max_precision;
}


/**
\brief Does not a thing, because cannot.

Cannot change precision of fixed precision hardware doubles.  This function is provided for template lookup, and interface completeness.

\param times Some times \f$\in \mathbb{C}\f$.
\param samples Some space samples \f$\in \mathbb{C}^n\f$.

\return The precision, which is now uniform.
*/
inline
unsigned EnsureAtUniformPrecision(TimeCont<dbl> & times, SampCont<dbl> & derivatives)
{
	return DoublePrecision();
}


/**
\brief Changes precision of mpfr to highest needed precision for the samples.

\param times Some times \f$\in \mathbb{C}\f$.
\param samples Some space samples \f$\in \mathbb{C}^n\f$.

\return The precision, which is now uniform.
*/
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples)
{
	auto def_prec = DefaultPrecision();
	if (std::any_of(begin(times),end(times),[=](auto const& p){return Precision(p)!=def_prec;}) 
		||
		std::any_of(begin(samples),end(samples),[=](auto const& p){return Precision(p)!=def_prec;}))
	{
		auto max_precision = max(MaxPrecision(samples), MaxPrecision(times));

		DefaultPrecision(max_precision);
		SetPrecision(times, max_precision);
		SetPrecision(samples, max_precision);
		return max_precision;
	}
	return def_prec;
	
}


/**
\brief Changes precision of mpfr to highest needed precision for the samples.

This function does NOT do any refinement, it merely changes the precision of default, and of the input objects.

\param times The times of some space samples
\param samples Some spatial samples, probably obtained from an Endgame
\param derivatives The derivatives of the space samples, at the time samples.  Again, probably obtained from Endgame.

\return The precision changed to.
*/
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples, SampCont<mpfr> & derivatives)
{
	auto def_prec = DefaultPrecision();
	if (std::any_of(begin(samples),end(samples),[=](auto const& p){return Precision(p)!=def_prec;}) 
	    || 
	    std::any_of(begin(derivatives),end(derivatives),[=](auto const& p){return Precision(p)!=def_prec;}))
	{
		auto max_precision = max(MaxPrecision(samples),MaxPrecision(times),MaxPrecision(derivatives));

		DefaultPrecision(max_precision);
		
		SetPrecision(times, max_precision);
		SetPrecision(samples, max_precision);
		SetPrecision(derivatives, max_precision);
		return max_precision;
	}
	return def_prec;
}



}}} // re: namespaces
