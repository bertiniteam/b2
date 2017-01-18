//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/trace.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/trace.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/trace.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake, Notre Dame


/**
\file bertini2/nag_algorithms/trace.hpp 

\brief Provides methods for computing traces.  

*/

#pragma once

#include "bertini2/nag_datatypes/witness_set.hpp"



namespace bertini {

	namespace algorithms {


		template <typename ComplexT>
		ComplexT Trace(nag_datatype::WitnessSet<ComplexT> const& w)
		{
			return ComplexT(0);
		}

	} // algorithm

} // bertini