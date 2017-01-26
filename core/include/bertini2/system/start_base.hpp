//This file is part of Bertini 2.
//
//start_base.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_base.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_base.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame

/**
\file b2/core/include/bertini2/system/start_base.hpp 

\brief Defines generic start system type.
*/


#pragma once

#include "bertini2/system/system.hpp"
#include "bertini2/limbo.hpp"


namespace bertini 
{
	namespace start_system{

		/**
		\brief Abstract base class for other start systems.

		Abstract base class for other start systems.  Start systems are special types of systems, to which we know solutions.  We also know how to construct various types of start systems from arbitrary polynomial systems.

		This class provides the empty virtual declarations for necessary override functions for specific start systems, including NumStartPoints (provides an upper bound on the number of solutions to the target system), and the private functions GenerateStartPoint(index), in double and multiple precision.  These two Generate functions are called by the templated non-overridden function StartPoint(index), which calls the appropriate one based on template type.
		*/
		class StartSystem : public System
		{

		public:
			

			virtual unsigned long long NumStartPoints() const = 0;

			template<typename T>
			Vec<T> StartPoint(unsigned long long index) const
			{
				return GenerateStartPoint(T(),index);
			}

			virtual ~StartSystem() = default;
			
		private:
			virtual Vec<dbl> GenerateStartPoint(dbl,unsigned long long index) const = 0;
			virtual Vec<mpfr> GenerateStartPoint(mpfr,unsigned long long index) const = 0;

			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<System>(*this);
			}

		};

	}
}



