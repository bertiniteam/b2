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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame

/**
\file bertini2/system/start_base.hpp 

\brief Defines generic start system type.
*/


#pragma once

#include "bertini2/system/system.hpp"

#include <functional>
#include <memory>


namespace bertini
{
	namespace start_system{

		// forward declare so the StartSystemFactory alias / MakeStartFactory below can be defined
		// here (they only need the StartSystem base, declared just below).
		class StartSystem;

		/**
		\brief A factory that builds a start system from a (prepared) target system.

		The zero-dim algorithm holds its start system polymorphically through the StartSystem base;
		the construction site (which knows the concrete start-system type) hands ZeroDimSolver one of
		these factories, which it calls on the homogenized/patched target to mint the start system.
		Adding a new start system therefore costs no new solver instantiation -- just a factory.
		*/
		template<typename SystemType>
		using StartSystemFactory =
			std::function<std::shared_ptr<StartSystem>(SystemType const&)>;

		/**
		\brief Build a StartSystemFactory that mints a concrete start system from the target.

		The single place a concrete start-system type is named when wiring up a ZeroDimSolver.
		`MakeStartFactory<start_system::TotalDegree>()` yields a factory that `make_shared`s a
		TotalDegree from the (prepared) target.  Used by the blackbox switch ladder, the python
		bindings, ZeroDimSolver's default, and the C++ tests.  This template is only instantiated
		where the concrete StartType is complete, so start_base.hpp need not include the concrete
		start systems.
		*/
		template<typename StartType, typename SystemType = bertini::System>
		StartSystemFactory<SystemType> MakeStartFactory()
		{
			return [](SystemType const& target) -> std::shared_ptr<StartSystem> {
				return std::make_shared<StartType>(target);
			};
		}

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
			virtual Vec<complex_dbl> GenerateStartPoint(complex_dbl,unsigned long long index) const = 0;
			virtual Vec<complex_mp> GenerateStartPoint(complex_mp,unsigned long long index) const = 0;

			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned /*version*/) {
				ar & boost::serialization::base_object<System>(*this);
			}

		};

	}
}



