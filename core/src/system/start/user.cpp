//This file is part of Bertini 2.
//
//mhom.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mhom.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mhom.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/start/user.hpp"


BOOST_CLASS_EXPORT(bertini::start_system::User);


namespace bertini {

	namespace start_system {

		// constructor for User start system, from any other *suitable* system.
		User::User(System const& s, SampCont<complex_dbl> const& solns) : user_system_(s), solns_in_dbl_(true)
		{
			std::get<SampCont<complex_dbl>>(solns_) = solns;
		}

		User::User(System const& s, SampCont<complex_mp> const& solns) : user_system_(s), solns_in_dbl_(false)
		{
			std::get<SampCont<complex_mp>>(solns_) = solns;
		}
				
		
		unsigned long long User::NumStartPoints() const
		{
			if (solns_in_dbl_)
				return std::get<SampCont<complex_dbl>>(solns_).size();
			else
				return std::get<SampCont<complex_mp>>(solns_).size();
		}


		
		Vec<complex_dbl> User::GenerateStartPoint(complex_dbl,unsigned long long index) const
		{
			if (solns_in_dbl_)
				return std::get<SampCont<complex_dbl>>(solns_)[index];
			else
			{
				const auto& r = std::get<SampCont<complex_mp>>(solns_)[index];
				Vec<complex_dbl> pt(r.size());
				for (unsigned ii=0; ii<r.size(); ++ii)
					pt(ii) = complex_dbl(r(ii));

				return pt;
			}
		}


		Vec<complex_mp> User::GenerateStartPoint(complex_mp,unsigned long long index) const
		{
			if (solns_in_dbl_)
			{
				const auto& r = std::get<SampCont<complex_dbl>>(solns_)[index];
				Vec<complex_mp> pt(r.size());
				for (unsigned ii=0; ii<r.size(); ++ii)
					pt(ii) = static_cast<complex_mp>(r(ii));

				return pt;
			}
			else
			{
				return std::get<SampCont<complex_mp>>(solns_)[index];
			}
		}

	} // namespace start_system
} //namespace bertini
