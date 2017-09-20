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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

#include "bertini2/system/start/user.hpp"


BOOST_CLASS_EXPORT(bertini::start_system::User);


namespace bertini {

	namespace start_system {

		// constructor for User start system, from any other *suitable* system.
		User::User(System const& s, SampCont<dbl> const& solns) : user_system_(s), solns_in_dbl_(true)
		{
			std::get<SampCont<dbl>>(solns_) = solns;
		}

		User::User(System const& s, SampCont<mpfr> const& solns) : user_system_(s), solns_in_dbl_(false)
		{
			std::get<SampCont<mpfr>>(solns_) = solns;
		}
				
		
		unsigned long long User::NumStartPoints() const
		{
			if (solns_in_dbl_)
				return std::get<SampCont<dbl>>(solns_).size();
			else
				return std::get<SampCont<mpfr>>(solns_).size();
		}


		
		Vec<dbl> User::GenerateStartPoint(dbl,unsigned long long index) const
		{
			if (solns_in_dbl_)
				return std::get<SampCont<dbl>>(solns_)[index];
			else
			{
				const auto& r = std::get<SampCont<mpfr>>(solns_)[index];
				Vec<dbl> pt(r.size());
				for (unsigned ii=0; ii<r.size(); ++ii)
					pt(ii) = dbl(r(ii));

				return pt;
			}
		}


		Vec<mpfr> User::GenerateStartPoint(mpfr,unsigned long long index) const
		{
			if (solns_in_dbl_)
			{
				const auto& r = std::get<SampCont<dbl>>(solns_)[index];
				Vec<mpfr> pt(r.size());
				for (unsigned ii=0; ii<r.size(); ++ii)
					pt(ii) = static_cast<mpfr>(r(ii));

				return pt;
			}
			else
			{
				return std::get<SampCont<mpfr>>(solns_)[index];
			}
		}

	} // namespace start_system
} //namespace bertini
