//This file is part of Bertini 2.
//
//python/zero_dim_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/zero_dim_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/zero_dim_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2023 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin-Eau Claire
//  2023
//
//
//  python/zero_dim_export.hpp:  Header file for exposing the zero dim solve algorithm to Python


#ifndef BERTINI2_PYBERTINI_NAG_ALGORITHMS_ZERO_DIM
#define BERTINI2_PYBERTINI_NAG_ALGORITHMS_ZERO_DIM

#pragma once

#include "python_common.hpp"

#include <bertini2/endgames.hpp>
#include <bertini2/nag_algorithms/zero_dim_solve.hpp>





namespace bertini{
	namespace python{
		
		using namespace bertini;




void ExportZeroDimAlgorithms();





template<typename AlgoT>
class ZDVisitor: public def_visitor<ZDVisitor<AlgoT> >
{


	friend class def_visitor_access;
		
	public:
		template<class PyClass>
		void visit(PyClass& cl) const;
		
	private:



};






}} // namespaces




#endif // the include guards
