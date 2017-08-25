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
		User::User(System const& s)
		{

			if (s.NumHomVariableGroups() > 0)
				throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");


			if (s.NumTotalFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct total degree start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct total degree start system, but target system has path varible declared already");			

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");
		}

		
		User& User::operator*=(Nd const& n)
		{
			*this *= n;
			return *this;
		}
		
		

		unsigned long long User::NumStartPoints() const
		{
			return 0;
		}


		
		Vec<dbl> User::GenerateStartPoint(dbl,unsigned long long index) const
		{
			Vec<dbl> start_point(NumVariables());

			return start_point;
		}


		Vec<mpfr> User::GenerateStartPoint(mpfr,unsigned long long index) const
		{
			Vec<mpfr> start_point(NumVariables());

			return start_point;
		}

		inline
		User operator*(User td, std::shared_ptr<node::Node> const& n)
		{
			td *= n;
			return td;
		}

	} // namespace start_system
} //namespace bertini
