// python/nid_datatypes_export.cpp — NID result data type registrations (WitnessSet, NIDResult).
// These are pure data types; no algorithm class_<> instantiations needed.

#include "numerical_irreducible_decomposition_export.hpp"
#include <boost/python/copy_const_reference.hpp>

namespace bertini{
	namespace python{

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

}} // namespaces
