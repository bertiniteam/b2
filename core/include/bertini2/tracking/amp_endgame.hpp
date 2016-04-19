//This file is part of Bertini 2.
//
//amp_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University



#pragma once

/**
\file base_endgame.hpp

\brief Contains parent class, Endgame, the parent class for all endgames.
*/

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/limbo.hpp"
#include "bertini2/tracking/tracking_config.hpp"
#include "bertini2/tracking/amp_tracker.hpp"

namespace bertini{ namespace tracking { namespace endgame {
			

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
unsigned EnsureAtUniformPrecision(TimeCont<dbl> & times, SampCont<dbl> & derivatives)
{
	return DoublePrecision();
}


//changes precision of mpfr to highest needed precision for the samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples)
{
	using std::max;
	unsigned max_precision = max(
	                             mpfr_float::default_precision(),
	                             MaxPrecision(samples)
	                             ); 

	if (max_precision != mpfr_float::default_precision())
		BOOST_LOG_TRIVIAL(severity_level::trace) << "EG changing default precision from " << mpfr_float::default_precision() << " to " << max_precision << std::endl;
	
	mpfr_float::default_precision(max_precision);

	SetPrecision(times, max_precision);
	SetPrecision(samples, max_precision);

	return max_precision;
}


//changes precision of mpfr to highest needed precision for the samples.
inline
unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples, SampCont<mpfr> & derivatives)
{
	unsigned max_precision = max(
	                             mpfr_float::default_precision(),
	                             MaxPrecision(samples),
	                             MaxPrecision(derivatives)
	                             ); 

	if (max_precision != mpfr_float::default_precision())
		BOOST_LOG_TRIVIAL(severity_level::trace) << "EG changing default precision from " << mpfr_float::default_precision() << " to " << max_precision << std::endl;
	
	mpfr_float::default_precision(max_precision);

	SetPrecision(times, max_precision);
	SetPrecision(samples, max_precision);
	SetPrecision(derivatives, max_precision);

	return max_precision;
}







/**
\brief Abstract class, defining a policy for Adaptive precision.

It is expected that this will be derived from for the various AMP instantiations of the endgames.
*/
class AMPEndgamePolicyBase
{
protected:
	
	// required virtual implementations
	virtual void ChangeTemporariesPrecisionImpl(unsigned new_precision) const = 0;
	virtual void MultipleToMultipleImpl(unsigned new_precision) const = 0;
	virtual void DoubleToMultipleImpl(unsigned new_precision) const = 0;
	virtual void MultipleToDoubleImpl() const = 0;



	mutable unsigned precision_;
	unsigned initial_precision_;
	bool preserve_precision_ = false; ///< Whether the endgame should change back to the initial precision after running 

	auto Precision() const
	{ return precision_; }



	/**
	\brief Change precision of all temporary internal state variables.

	This excludes those which cannot be re-written without copying.

	\brief new_precision The new precision to adjust to.
	*/
	void ChangeTemporariesPrecision(unsigned new_precision) const
	{
		ChangeTemporariesPrecisionImpl(new_precision);
	}

	/**
	\brief Converts from multiple to different precision multiple precision

	Changes the precision of the internal temporaries to desired precision

	\param new_precision The new precision.
	*/
	void MultipleToMultiple(unsigned new_precision) const
	{
		mpfr_float::default_precision(new_precision);
		precision_ = new_precision;

		ChangeTemporariesPrecision(new_precision);

		MultipleToMultipleImpl(new_precision);
	}

	/**
	\brief Converts from multiple to different precision multiple precision

	Changes the precision of the internal temporaries to desired precision

	\param new_precision The new precision.
	*/
	void DoubleToMultiple(unsigned new_precision) const
	{
		mpfr_float::default_precision(new_precision);
		precision_ = new_precision;
		
		ChangeTemporariesPrecision(new_precision);

		DoubleToMultipleImpl(new_precision);
	}

	/**
	\brief Converts from multiple to different precision multiple precision

	Changes the precision of the internal temporaries to desired precision

	\param new_precision The new precision.
	*/
	void MultipleToDouble() const
	{
		mpfr_float::default_precision(DoublePrecision());
		precision_ = DoublePrecision();

		ChangeTemporariesPrecision(DoublePrecision());

		MultipleToDoubleImpl();
	}

	bool PrecisionSanityCheck() const
	{
		return true;
	}

public:	

	/**
	Change precision of endgame object to next_precision.  Converts the internal temporaries, and adjusts precision of system. 

	\param new_precision The precision to change to.
	\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
	*/
	SuccessCode ChangePrecision(unsigned new_precision) const
	{
		if (new_precision==precision_) // no op
			return SuccessCode::Success;

		if (new_precision==DoublePrecision() && precision_>DoublePrecision())
			MultipleToDouble();
		else if(new_precision > DoublePrecision() && precision_ == DoublePrecision())
			DoubleToMultiple(new_precision);
		else
			MultipleToMultiple(new_precision);


		#ifndef BERTINI_DISABLE_ASSERTS
		assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
		#endif

		return SuccessCode::Success;
	}

	AMPEndgamePolicyBase() : precision_(mpfr_float::default_precision())
	{}

	virtual ~AMPEndgamePolicyBase() = default;
}; // re: class AMPEndgamePolicyBase
			

} } } // end namespaces 
				


