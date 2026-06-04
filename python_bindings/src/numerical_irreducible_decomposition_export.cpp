//This file is part of Bertini 2.
//
//python/numerical_irreducible_decomposition_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/numerical_irreducible_decomposition_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/numerical_irreducible_decomposition_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/numerical_irreducible_decomposition_export.cpp:  source file for exposing the numerical irreducible decomposition algorithm to python.

#include "numerical_irreducible_decomposition_export.hpp"
#include "configured_visitor.hpp"

#include <boost/python/copy_const_reference.hpp>



namespace bertini{
	namespace python{

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




void ExportNID(){
	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".nag_algorithms");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("nag_algorithms") = new_submodule;


	scope new_submodule_scope = new_submodule;

	// the shared config structs (Tolerances, Regeneration, Sharpening, PostProcessing)
	// are registered by ExportZeroDim, which runs before this.  We reuse them; the
	// ConfiguredVisitor finds them via the boost.python registry.
	ExportNIDDataTypes();
	ExportNIDAlgorithms();

}




// a helper for exposing a single concrete witness-set type.
template<typename NumT>
void ExportWitnessSet(std::string const& class_name){

	using WS = nag_datatype::WitnessSet<NumT>;

	class_<WS>(class_name.c_str(), init<>())
	.def("degree", &WS::Degree, "the degree of the component, i.e. the number of witness points")
	.def("dimension", &WS::Dimension, "the dimension of the component")
	.def("is_consistent", &WS::IsConsistent, "whether the slice dimension matches the system's underdeterminedness")
	.def("get_point", &WS::GetPoint, return_value_policy<copy_const_reference>(), "get the i-th witness point")
	.def("get_system", &WS::GetSystem, return_internal_reference<>(), "the system this witness set is for")
	;
}



// a helper for exposing a single concrete NID result type.
template<typename NumT>
void ExportNIDResult(std::string const& class_name){

	using R = nag_datatype::NumericalIrreducibleDecomposition<NumT>;

	class_<R>(class_name.c_str(), init<>())
	.def("nonempty_codimensions", &R::NonEmptyCodimensions, "the distinct codimensions which contain at least one component")
	.def("num_witness_sets", &R::NumWitnessSets, "the number of stored witness sets")
	.def("get_witness_set", &R::GetWitnessSet, return_internal_reference<>(), "get the i-th stored witness set")
	;
}



void ExportNIDDataTypes(){

	ExportWitnessSet<dbl_complex>("WitnessSetDoublePrecision");
	ExportWitnessSet<mpfr_complex>("WitnessSetMultiplePrecision");

	ExportNIDResult<dbl_complex>("NumericalIrreducibleDecompositionDoublePrecision");
	ExportNIDResult<mpfr_complex>("NumericalIrreducibleDecompositionMultiplePrecision");
}




// a helper function used immediately below.  There is no declaration in the header file...
template<typename TrackerT, typename EndgameT>
void ExportNIDSpecific(std::string const& class_name){

	using NIDT = algorithm::NumericalIrreducibleDecomposition<TrackerT, EndgameT, bertini::System>;

	class_<NIDT, std::shared_ptr<NIDT> >(class_name.c_str(), init<bertini::System>())
	.def(NIDVisitor<NIDT>())
	;
}




void ExportNIDAlgorithms(){

	{
		using TrackerT = bertini::tracking::DoublePrecisionTracker;
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesDoublePrecision");
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyDoublePrecision");
	}

	{
		using TrackerT = bertini::tracking::MultiplePrecisionTracker;
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesFixedMultiplePrecision");
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyFixedMultiplePrecision");
	}

	{
		using TrackerT = bertini::tracking::AMPTracker;
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesAdaptivePrecision");
		ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyAdaptivePrecision");
	}

}



}} // namespaces
