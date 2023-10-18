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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file bertini2/system/start/user.hpp

\brief Defines the User-defined start system type.

the user provided start system is merely their system, evaluated at the start time.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/system/start/utility.hpp"
#include "bertini2/common/config.hpp"

// these next two blobs solve a problem of private-ness.
// see https://stackoverflow.com/questions/4112075/trouble-overriding-save-construct-data-when-serializing-a-pointer-to-a-class-wit

// forward declare the User start system
namespace bertini{ namespace start_system{
	class User;
}}

// forward declare the function, so we can friend it below.
namespace boost { namespace serialization {
template<class Archive>
inline void save_construct_data(Archive & ar, const bertini::start_system::User * t, const unsigned int file_version);
}}
// end nonsense for friends.  so lonely, but c++ friends don't solve the irl problem at all.  

namespace bertini
{
	namespace start_system{


		/**
		\brief The user-provided start system for Numerical Algebraic Geometry


		*/
		class User : public StartSystem
		{
		public:
			User() = delete; // deleted because requires a reference to a System to construct
			virtual ~User() = default;

			/**
			 Constructor for making a user-provided start system from another.
			*/
			User(System const& s, SampCont<dbl> const& solns);
			User(System const& s, SampCont<mpfr_complex> const& solns);



			/**
			Get the number of start points for this start system.
			*/
			unsigned long long NumStartPoints() const override;

			User& operator*=(Nd const& n) = delete;

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
			Vec<mpfr_complex> GenerateStartPoint(mpfr_complex,unsigned long long index) const override;


			friend class boost::serialization::access;
			template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const User * t, const unsigned int file_version);

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<StartSystem>(*this);
			}

			const bertini::System& user_system_;
			std::tuple<SampCont<dbl>, SampCont<mpfr_complex>> solns_;
			bool solns_in_dbl_;
		};
	}
}

namespace boost { namespace serialization {
template<class Archive>
inline void save_construct_data(
    Archive & ar, const bertini::start_system::User * t, const unsigned int file_version
){
    // save data required to construct instance
    ar << t->user_system_;
    ar << t->solns_in_dbl_;

    ar << std::get<0>(t->solns_);
	ar << std::get<1>(t->solns_);

}

template<class Archive>
inline void load_construct_data(
    Archive & ar, bertini::start_system::User * t, const unsigned int file_version
){
    // retrieve data from archive required to construct new instance
    bertini::System sys;
    ar >> sys;

    bool solns_in_dbl;
    ar >> solns_in_dbl;

    if (solns_in_dbl)
    {
    	bertini::SampCont<bertini::dbl> solns;
    	ar >> solns;
    	::new(t)bertini::start_system::User(sys, solns);
    }
	else
	{
		bertini::SampCont<bertini::mpfr_complex> solns;
		ar >> solns;
		::new(t)bertini::start_system::User(sys, solns);
	}
}
}}
