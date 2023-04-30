//This file is part of Bertini 2.
//
//prec_base.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//prec_base.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with prec_base.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Tim Hodges, Colorado State University



#pragma once


/**
\file include/bertini2/endgames/prec_base.hpp

\brief Contains a parent class, EndgamePrecPolicyBase (an observable), from which the fixed double, fixed multiple, or adaptive precision endgames are derived.
*/


namespace bertini{ namespace endgame {


/**
\brief A common base type for various precision types, fixed and adaptive.  The purpose of this is to maintain a uniform interface to the tracker that's being used, across endgame types.
*/
template <typename TrackerT>
class EndgamePrecPolicyBase : public virtual Observable
{
public:

	using TrackerType = TrackerT;

	explicit
	EndgamePrecPolicyBase(TrackerT const& new_tracker) : tracker_(std::ref(new_tracker))
	{}

	virtual ~EndgamePrecPolicyBase() = default;
	
	/**
	Tell the endgame to use the given tracker.  Takes a reference.  

	\note Ensure the tracker you are using doesn not go out of scope!
	*/
	inline
	void SetTracker(TrackerT const& new_tracker)
	{
		tracker_ = std::ref(new_tracker); // rebind the reference
	}


	/**
	\brief Getter for the tracker used inside an instance of the endgame. 
	*/
	inline
	const TrackerT& GetTracker() const
	{
		return tracker_.get();
	}

	/**
	\brief Get the system being tracked on, which is referred to by the tracker.
	*/
	inline
	const System& GetSystem() const 
	{ return GetTracker().GetSystem();}

	void ChangePrecision(unsigned p)
	{
		tracker_.ChangePrecision(p);
	}
	
private:

	/**
	\brief A tracker that must be passed into the endgame through a constructor. This tracker is what will be used to track to all time values during the endgame. 
	*/
	std::reference_wrapper<const TrackerT> tracker_;

}; //EndgamePrecPolicyBase



} } // end namespaces 
				


