//This file is part of Bertini 2.
//
//nag_algorithms/common/algorithm_base.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nag_algorithms/common/algorithm_base.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nag_algorithms/common/algorithm_base.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file nag_algorithms/common/algorithm_base.hpp

\brief provides a common base type for algorithms to inherit from, for the abstraction layer above for calling them from blackbox-type programs
*/

#pragma once

namespace bertini{

namespace algorithm {

	struct AnyAlgorithm
	{

		virtual 
		void Run() = 0;

		virtual ~AnyAlgorithm() = default;
	};
}
}

