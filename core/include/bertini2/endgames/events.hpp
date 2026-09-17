//This file is part of Bertini 2.
//
//events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


/**
\file endgames/events.hpp

\brief Contains the endgames/events base types
*/

#pragma once

#include "bertini2/detail/events.hpp"

namespace bertini {

	namespace endgame{


	
	/**
	\brief Generic event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameEvent,ConstEvent);

	/**
	\brief Generic failure event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameFailure,ConstEvent);

	/**
	\brief Generic success event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameSuccess,ConstEvent);

	/**
	\brief MinTrackTime reached during endgame
	*/
	ADD_BERTINI_EVENT_TYPE(MinTrackTimeReached,EndgameFailure);

	/**
	\brief Refining a sample failed for some reason
	*/
	ADD_BERTINI_EVENT_TYPE(RefiningFailed,EndgameFailure);

	/**
	\brief Cycle number computed was too high
	*/
	ADD_BERTINI_EVENT_TYPE(CycleNumTooHigh,EndgameFailure);

	/**
	\brief Security max norm reached. 
	*/
	ADD_BERTINI_EVENT_TYPE(SecurityMaxNormReached,EndgameFailure);

	/**
	\brief Started running the endgame
	*/
	ADD_BERTINI_EVENT_TYPE(Initializing,EndgameEvent);

	/**
	\brief The adaptive endgame needed a higher precision before it had its first
	approximation, and recomputes its sample window at the new precision.

	The samples it announced at the lower precision are superseded -- the endgame discards
	them and tracks the window again from the boundary point -- so an observer keeping the
	sequence of samples should drop them too.  (A precision increase LATER in the run, once
	approximations exist, keeps the samples and widens them in place; that one announces only
	PrecisionChanged.)  Emitted right after the matching PrecisionChanged.
	*/
	ADD_BERTINI_EVENT_TYPE(SamplesRecomputedAtHigherPrecision,EndgameEvent);


	/**
	\brief Time advancing
	*/
	ADD_BERTINI_EVENT_TYPE(TimeAdvanced,EndgameEvent);

	/**
	\brief A new sample point on the path was computed, adding one point to the sequence of the
	endgame's approach to the target time.

	Emitted by BOTH endgame flavors, wherever a new path sample becomes final: the
	power series endgame's time advance (after the new sample is refined), and the
	Cauchy endgame's rotation onto its power-series window (as tracked -- Cauchy does
	not refine at that site).  Together the emissions form a time-indexed sequence of
	points approaching the root, which is what lets a consumer watch a derived
	quantity -- the singular values of a Jacobian, say -- behave as a function of
	distance to the root, rather than thresholding it at a single point.

	The sample and its time ride ON the event, as with CircleAdvanced, so an observer
	never has to know which flavor emitted it: the two flavors keep their samples in
	differently-named containers (GetSamples() versus GetPSEGSamples()).

	Distinct from TimeAdvanced, which announces only that time moved and carries no
	payload, and from SampleRefined, which fires once per EXISTING sample when the
	whole window is re-refined rather than once per new sample.
	*/
	template<class ObservedT>
	class ComputedSamplePoint : public EndgameEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		using ComplexT = typename ObservedT::BaseComplexT;  ///< The complex number type.

		/**
		\brief The constructor for a ComputedSamplePoint event.

		\param obs The observable emitting the event.
		\param new_point The newly computed sample point on the path.
		\param new_time The time value at which the sample was computed.
		*/
		ComputedSamplePoint(const ObservedT & obs,
		                    Vec<ComplexT> const& new_point,
		                    ComplexT const& new_time) : EndgameEvent<ObservedT>(obs),
		                                                new_point_(new_point),
		                                                new_time_(new_time)
		{}

		virtual ~ComputedSamplePoint() = default;
		ComputedSamplePoint() = delete;

		/// \return The newly computed sample point on the path.
		const auto& NewSample() const {return new_point_;}

		/// \return The time value at which the sample was computed.
		const auto& NewTime() const {return new_time_;}

	private:
		const Vec<ComplexT>& new_point_;
		const ComplexT& new_time_;
	};


	/**
	\brief Advanced around the circle around target time
	*/
	template<class ObservedT>
	class CircleAdvanced : public EndgameEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		using ComplexT = typename ObservedT::BaseComplexT;  ///< The complex number type.
		/**
		\brief The constructor for a CircleAdvanced Event.

		\param obs The observable emitting the event.
		\param new_point The space point reached after advancing around the circle.
		\param new_time The time value reached after advancing around the circle.
		*/
		CircleAdvanced(const ObservedT & obs,
		               Vec<ComplexT> const& new_point,
		               ComplexT const& new_time) : EndgameEvent<ObservedT>(obs),
													new_point_(new_point),
													new_time_(new_time)
		{}


		virtual ~CircleAdvanced() = default;
		CircleAdvanced() = delete;
		
		/// \return The space point reached after advancing around the circle.
		const auto& NewSample() const {return new_point_;}

		/// \return The time value reached after advancing around the circle.
		const auto& NewTime() const {return new_time_;}

	private:
		const Vec<ComplexT>& new_point_;
		const ComplexT& new_time_;
	};


	/**
	\brief Walked a complete loop around the target time.
	*/
	ADD_BERTINI_EVENT_TYPE(ClosedLoop,EndgameEvent);

	/**
	\brief Approximated a root at target time.
	*/
	ADD_BERTINI_EVENT_TYPE(ApproximatedRoot,EndgameEvent);

	/**
	\brief Converged -- endgame is done!
	*/
	ADD_BERTINI_EVENT_TYPE(Converged,EndgameSuccess);

	/**
	\brief Refined a sample
	*/
	ADD_BERTINI_EVENT_TYPE(SampleRefined,EndgameEvent);

	/**
	\brief Made it into the EG operating zone, or so we believe
	*/
	ADD_BERTINI_EVENT_TYPE(InEGOperatingZone,EndgameEvent);
	

	/// \brief Event emitted when the endgame's working precision changes.
	template<class ObservedT>
	class PrecisionChanged : public EndgameEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a PrecisionChanged Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		PrecisionChanged(const ObservedT & obs, 
		             unsigned previous, unsigned next) : EndgameEvent<ObservedT>(obs),
													prev_(previous),
													next_(next)
		{}


		virtual ~PrecisionChanged() = default;
		PrecisionChanged() = delete;
		
		/**
		\brief Get the previous precision.
		*/
		auto Previous() const {return prev_;}

		/**
		\brief Get the next precision, what it changed to.
		*/
		auto Next() const {return next_;}
	private:
		const unsigned prev_, next_;
	};


	}// re: namespace endgames
}// re: namespace bertini

