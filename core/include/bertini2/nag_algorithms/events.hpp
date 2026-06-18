//This file is part of Bertini 2.
//
//nag_algorithms/events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nag_algorithms/events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nag_algorithms/events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file nag_algorithms/events.hpp

\brief Events emitted by the higher-level algorithms (e.g. ZeroDim) during a solve.

These are emitted on a non-templated emitter base (e.g. AnyZeroDim) so that a single
observer type can watch every templated instantiation of the algorithm.  An observer
recovers the concrete algorithm (and its full interface) from the event's Get(): in
C++ via dynamic_cast, in Python automatically via RTTI to the most-derived registered
type.
*/

#pragma once

#include "bertini2/detail/events.hpp"

#include <cstddef>

namespace bertini {

	namespace algorithm {

		/**
		\brief Generic event for a numerical algebraic geometry algorithm.
		*/
		ADD_BERTINI_EVENT_TYPE(AlgorithmEvent, ConstEvent);

		/**
		\brief The algorithm has begun (e.g. the whole zero-dim solve is starting).
		*/
		ADD_BERTINI_EVENT_TYPE(AlgorithmStarted, AlgorithmEvent);

		/**
		\brief The algorithm has finished.
		*/
		ADD_BERTINI_EVENT_TYPE(AlgorithmComplete, AlgorithmEvent);


		/**
		\brief The algorithm is beginning work on one solution path, carrying its index.

		Between this and the matching PathComplete, the algorithm runs everything for
		that path -- the main homotopy track and any endgame sub-tracks -- so an
		observer can group them together as one complete solution path.
		*/
		template<class ObservedT>
		class PathStarted : public AlgorithmEvent<ObservedT>
		{ BOOST_TYPE_INDEX_REGISTER_CLASS
		public:
			using HeldT = typename AlgorithmEvent<ObservedT>::HeldT;
			PathStarted(HeldT obs, std::size_t path_index)
				: AlgorithmEvent<ObservedT>(obs), path_index_(path_index)
			{}

			std::size_t PathIndex() const { return path_index_; }

			virtual ~PathStarted() = default;
			PathStarted() = delete;
		private:
			std::size_t path_index_;
		};

		/**
		\brief The algorithm has finished work on one solution path, carrying its index.
		*/
		template<class ObservedT>
		class PathComplete : public AlgorithmEvent<ObservedT>
		{ BOOST_TYPE_INDEX_REGISTER_CLASS
		public:
			using HeldT = typename AlgorithmEvent<ObservedT>::HeldT;
			PathComplete(HeldT obs, std::size_t path_index)
				: AlgorithmEvent<ObservedT>(obs), path_index_(path_index)
			{}

			std::size_t PathIndex() const { return path_index_; }

			virtual ~PathComplete() = default;
			PathComplete() = delete;
		private:
			std::size_t path_index_;
		};

	} // namespace algorithm

} // namespace bertini
