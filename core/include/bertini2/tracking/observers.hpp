//This file is part of Bertini 2.0.
//
//tracking/observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking/observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  tracking/observers.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file tracking/observers.hpp

\brief Contains the tracking/observers base types
*/

#pragma once

#include "bertini2/tracking/events.hpp"
#include "bertini2/tracking/base_tracker.hpp"
#include "bertini2/tracking/events.hpp"

namespace bertini {

	namespace tracking{

		template<class ObservedT>
		class PrecisionAccumulator : public Observer<ObservedT>
		{
			virtual void Observe(EventBase const& e) override
			{
				const TrackingEvent<ObservedT>* p = dynamic_cast<const TrackingEvent<ObservedT>*>(&e);
				if (p)
				{
					Visit(p->Get());
				}
				else
				{
					std::cout << "failed conversion to known event type...\n";
				}
			}


			virtual void Visit(ObservedT const& t) override
			{
				precisions_.push_back(t.CurrentPrecision());
			}

		public:
			const std::vector<unsigned>& Precisions() const
			{
				return precisions_;
			}

		private:
			std::vector<unsigned> precisions_;
		};

		/**
		Example usage:
		PathAccumulator<AMPTracker> path_accumulator;
		*/
		template<class ObservedT, template<class> class EventT = SuccessfulStep>
		class PathAccumulator : public Observer<ObservedT>
		{
			virtual void Observe(EventBase const& e) override
			{
				const EventT<ObservedT>* p = dynamic_cast<const EventT<ObservedT>*>(&e);
				if (p)
				{
					Visit(p->Get());
				}
				else
				{
					std::cout << "failed conversion to known event type...\n";
				}
			}


			virtual void Visit(ObservedT const& t) override
			{
				if (t.CurrentPrecision()==DoublePrecision())
				{
					Vec<mpfr> temp(t.NumVariables());
					for (unsigned ii=0; ii<t.NumVariables(); ++ii)
						temp(ii) = mpfr(t.template CurrentPoint<dbl>()(ii));
					
					path_.push_back(temp);
				}
				else
					path_.push_back(t.template CurrentPoint<mpfr>());
			}

		public:
			const std::vector<Vec<mpfr> >& Path() const
			{
				return path_;
			}

		private:
			std::vector<Vec<mpfr> > path_;
		};
	}
}// re: namespace bertini
