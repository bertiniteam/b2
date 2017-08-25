//This file is part of Bertini 2.
//
//bertini2/system/start/user.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/system/start/user.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/system/start/user.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file bertini2/system/start/user.hpp 

\brief Defines the User-defined start system type.

the user provided start system is merely their system, evaluated at the start time.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/limbo.hpp"


namespace bertini 
{
	namespace start_system{


		/**
		\brief The user-provided start system for Numerical Algebraic Geometry

		
		*/
		class User : public StartSystem
		{
		public:
			User() = default;
			virtual ~User() = default;
			
			/**
			 Constructor for making a user-provided start system from another.
			*/
			User(System const& s);


			

			/**
			Get the number of start points for this start system.
			*/
			unsigned long long NumStartPoints() const override;

			User& operator*=(Nd const& n);

			User& operator+=(System const& sys) = delete;
			
		private:

			/**
			Get the ith start point, in double precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<dbl> GenerateStartPoint(dbl,unsigned long long index) const override;

			/**
			Get the ith start point, in current default precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<mpfr> GenerateStartPoint(mpfr,unsigned long long index) const override;


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<StartSystem>(*this);
			}

		};
	}
}



