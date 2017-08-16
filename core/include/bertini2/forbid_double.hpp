//This file is part of Bertini 2.
//
//b2/core/include/forbid_double.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/core/include/forbid_double.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/core/include/forbid_double.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

#ifndef BERTINI_FORBID_MIXED_ARITHMETIC_WITH_DOUBLES_HPP
#define BERTINI_FORBID_MIXED_ARITHMETIC_WITH_DOUBLES_HPP

#pragma once

namespace boost{
	namespace multiprecision
	{
template <class Backend>
struct is_compatible_arithmetic_type<double, number<Backend> > : public mpl::false_ {};
	}
}


#endif
