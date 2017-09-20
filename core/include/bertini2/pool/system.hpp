//This file is part of Bertini 2.
//
//include/bertini2/pool/system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/pool/system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/pool/system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

// include/bertini2/pool/system.hpp:  provides the bertini::SystemPool class.


/**
\file include/bertini2/pool/system.hpp 
*/



#ifndef BERTINI_SYSTEM_POOL_HPP
#define BERTINI_SYSTEM_POOL_HPP

#include "bertini2/pool/pool.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {

	class SystemPool : public detail::Pool<System>
	{

	};

	template<typename NumT>
	class PointPool : public detail::Pool<Vec<NumT> >
	{

	};


} // re: namespace bertini

#endif
