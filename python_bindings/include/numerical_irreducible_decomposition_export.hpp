//This file is part of Bertini 2.
//
//python/numerical_irreducible_decomposition_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/numerical_irreducible_decomposition_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/numerical_irreducible_decomposition_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin-Eau Claire
//  2026
//
//
//  python/numerical_irreducible_decomposition_export.hpp:  Header file for exposing the numerical irreducible decomposition algorithm to Python


#ifndef BERTINI2_PYBERTINI_NAG_ALGORITHMS_NID
#define BERTINI2_PYBERTINI_NAG_ALGORITHMS_NID

#pragma once

#include "python_common.hpp"
#include "configured_visitor.hpp"

#include <bertini2/endgames.hpp>
#include <bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp>




namespace bertini{
	namespace python{

		using namespace bertini;




void ExportNID();

// sub-functions (defined in nid_{datatypes,double,mp,amp}_export.cpp)
void ExportNIDDataTypes();
void ExportNIDDouble();
void ExportNIDMP();
void ExportNIDAMP();


template<typename AlgoT>
class NIDVisitor: public def_visitor<NIDVisitor<AlgoT> >
{
	friend class ::boost::python::def_visitor_access;

	public:
		template<class PyClass>
		void visit(PyClass& cl) const;

	private:

		using MutableTrackerGetter = typename AlgoT::TrackerT& (AlgoT::*)();
		static MutableTrackerGetter GetTrackerMutable()
		{
			return &AlgoT::GetTracker;
		};

		using MutableEndgameGetter = typename AlgoT::EndgameT& (AlgoT::*)();
		static MutableEndgameGetter GetEndgameMutable()
		{
			return &AlgoT::GetEndgame;
		};
};


// Visitor body — template member function must be in header so each split TU can instantiate it.
template<typename AlgoT>
template<class PyClass>
void NIDVisitor<AlgoT>::visit(PyClass& cl) const
{
	cl
	.def(ConfiguredVisitor<AlgoT>())
	.def("solve", &AlgoT::Solve, "run the numerical irreducible decomposition with the currently stored settings (not yet implemented)")
	.def("get_tracker", GetTrackerMutable(), return_internal_reference<>(), "get a mutable reference to the Tracker being used")
	.def("get_endgame", GetEndgameMutable(), return_internal_reference<>(), "get a mutable reference to the Endgame being used")
	.def("decomposition", &AlgoT::GetDecomposition, return_internal_reference<>(), "get the most recently computed numerical irreducible decomposition")
	;
}


// Helper template — defined here so all split TUs can use it.
template<typename TrackerT, typename EndgameT>
void ExportNIDSpecific(std::string const& class_name){
	using NIDT = algorithm::NumericalIrreducibleDecomposition<TrackerT, EndgameT, bertini::System>;
	class_<NIDT, std::shared_ptr<NIDT> >(class_name.c_str(), init<bertini::System>())
	.def(NIDVisitor<NIDT>())
	;
}




}} // namespaces




#endif // the include guards
