//This file is part of Bertini 2.
//
//python/zero_dim_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/zero_dim_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/zero_dim_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/zero_dim_export.cpp:  source file for exposing the zero dim algorithm to python.

#include "zero_dim_export.hpp"



namespace bertini{
	namespace python{

template<typename AlgoT>
template<class PyClass>
void ZDVisitor<AlgoT>::visit(PyClass& cl) const
{
	cl
	.def("solve", &AlgoT::Solve, "run the zero dim algorithm with currently stored settings")
	;
}




void ExportZeroDimAlgorithms(){

	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".nag_algorithms");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("nag_algorithms") = new_submodule;
	

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "Algorithms for computing things, like point solutions to square systems (zerodim algorithm).";


	using TrackerT = bertini::tracking::DoublePrecisionTracker;

	using ZeroDimTotalDegree = algorithm::ZeroDim<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, start_system::TotalDegree>;

	class_<ZeroDimTotalDegree, std::shared_ptr<ZeroDimTotalDegree> >("ZeroDimDoublePrecisionTotalDegree", init<bertini::System>())
	.def(ZDVisitor<ZeroDimTotalDegree>())
	;
}




}} // namespaces