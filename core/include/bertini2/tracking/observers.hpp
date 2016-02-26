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

		// class TrackingObserver :
		// 	public Observer<TrackingObserver, Tracker>
		// {

		// };

		// public VisitorBase,
		class PathAccumulator : public Observer<Tracker<AMPTracker> >
		{
			virtual void Update(EventBase const& e) override
			{
				const TrackingEvent<Tracker<AMPTracker> >* p = dynamic_cast<const TrackingEvent<Tracker<AMPTracker> >*>(&e);
				if (p)
				{
					std::cout << "asdf\n";
					Visit(p->Get());
				}
				else
				{
					std::cout << "failed conversion to SuccessfulStep...\n";
				}
			}


			virtual void Visit(Tracker<AMPTracker> const& t) override
			{
				std::cout << "visiting tracker\n";
			}
		};

	}
}// re: namespace bertini
