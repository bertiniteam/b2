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
	.def("get_tracker", GetTrackerMutable(), return_internal_reference<>(), "get a mutable reference to the Tracker being used")
	.def("get_endgame", GetEndgameMutable(), return_internal_reference<>(), "get a mutable reference to the Endgame being used")
	.def("solutions", &AlgoT::FinalSolutions, return_internal_reference<>(), "get the solutions at the target time")
	.def("solution_metadata", &AlgoT::FinalSolutionMetadata, return_internal_reference<>(), "get the metadata for the solutions at the target time")
	.def("endgame_boundary_data", &AlgoT::EndgameBoundaryData, return_internal_reference<>(), "get the data for the state at the endgame boundary (when we switch from regular tracking to endgame tracking")

	;
}







void ExportZeroDim(){
	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".nag_algorithms");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("nag_algorithms") = new_submodule;
	

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "Algorithms for computing things, like point solutions to square systems (zerodim algorithm).";


	ExportZDConfigs();
	ExportZDAlgorithms();

}




// a helper function used immediately below.  There is no declaration in the header file...
template<typename TrackerT, typename EndgameT, typename SystemT, typename StartSystemT>
void ExportZeroDimSpecific(std::string const& class_name){

	using ZeroDimT = algorithm::ZeroDim<TrackerT, EndgameT, SystemT, StartSystemT>;

	class_<ZeroDimT, std::shared_ptr<ZeroDimT> >(class_name.c_str(), init<SystemT>())
	.def(ZDVisitor<ZeroDimT>())
	;
}





void ExportZDAlgorithms(){



	{
		using TrackerT = bertini::tracking::DoublePrecisionTracker;
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesDoublePrecisionTotalDegree");
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyDoublePrecisionTotalDegree");
	}

	{
		using TrackerT = bertini::tracking::MultiplePrecisionTracker;
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesFixedMultiplePrecisionTotalDegree");
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyFixedMultiplePrecisionTotalDegree");
	}

	{
		using TrackerT = bertini::tracking::AMPTracker;
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesAdaptivePrecisionTotalDegree");
		ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyAdaptivePrecisionTotalDegree");
	}

}



void ExportZDConfigs()
{
	using namespace bertini::algorithm;

	class_<TolerancesConfig>("TolerancesConfig", init<>())
	.def_readwrite("newton_before_endgame", &TolerancesConfig::newton_before_endgame)
	.def_readwrite("newton_during_endgame", &TolerancesConfig::newton_during_endgame)
	.def_readwrite("final_tolerance", &TolerancesConfig::final_tolerance)
	.def_readwrite("path_truncation_threshold", &TolerancesConfig::path_truncation_threshold)
	;


	class_<MidPathConfig>("MidPathConfig", init<>())
	.def_readwrite("same_point_tolerance", &MidPathConfig::same_point_tolerance)
	;


	class_<AutoRetrackConfig>("AutoRetrackConfig", init<>())
	.def_readwrite("midpath_decrease_tolerance_factor", &AutoRetrackConfig::midpath_decrease_tolerance_factor)
	;


	class_<SharpeningConfig>("SharpeningConfig", init<>())
	.def_readwrite("sharpendigits", &SharpeningConfig::sharpendigits)
	.def_readwrite("function_residual_tolerance", &SharpeningConfig::function_residual_tolerance)
	.def_readwrite("ratio_tolerance", &SharpeningConfig::ratio_tolerance)
	;

	class_<RegenerationConfig>("RegenerationConfig", init<>())
	.def_readwrite("remove_infinite_endpoints", &RegenerationConfig::remove_infinite_endpoints)
	.def_readwrite("higher_dimension_check", &RegenerationConfig::higher_dimension_check)
	.def_readwrite("start_level", &RegenerationConfig::start_level)
	.def_readwrite("newton_before_endgame", &RegenerationConfig::newton_before_endgame)
	.def_readwrite("newton_during_endgame", &RegenerationConfig::newton_during_endgame)
	.def_readwrite("final_tolerance", &RegenerationConfig::final_tolerance)
	;


	class_<PostProcessingConfig>("PostProcessingConfig", init<>())
	.def_readwrite("real_threshold", &PostProcessingConfig::real_threshold)
	.def_readwrite("endpoint_finite_threshold", &PostProcessingConfig::endpoint_finite_threshold)
	.def_readwrite("same_point_tolerance", &PostProcessingConfig::same_point_tolerance)
	;

	class_<ZeroDimConfig<dbl_complex>>("ZeroDimConfigDoublePrec", init<>())
	.def_readwrite("start_time", &ZeroDimConfig<dbl_complex>::start_time)
	.def_readwrite("target_time", &ZeroDimConfig<dbl_complex>::target_time)
	.def_readwrite("endgame_boundary", &ZeroDimConfig<dbl_complex>::endgame_boundary)
	;

	class_<ZeroDimConfig<mpfr_complex>>("ZeroDimConfigMultiprec", init<>())
	.def_readwrite("start_time", &ZeroDimConfig<mpfr_complex>::start_time)
	.def_readwrite("target_time", &ZeroDimConfig<mpfr_complex>::target_time)
	.def_readwrite("endgame_boundary", &ZeroDimConfig<mpfr_complex>::endgame_boundary)
	;

}



template<typename NumT>
void ExposeSolutionMetaData(std::string const& class_name){
	using namespace bertini::algorithm;

	using MDT = SolutionMetaData<NumT>;
	class_<MDT>(class_name.c_str(),init<>())
	.def_readwrite("path_index",&MDT::path_index)
	.def_readwrite("solution_index",&MDT::solution_index)

	.def_readwrite("precision_changed",&MDT::precision_changed)
	.def_readwrite("time_of_first_prec_increase",&MDT::time_of_first_prec_increase)
	.def_readwrite("max_precision_used",&MDT::max_precision_used)

	.def_readwrite("pre_endgame_success",&MDT::pre_endgame_success)

	.def_readwrite("condition_number",&MDT::condition_number)
	.def_readwrite("newton_residual",&MDT::newton_residual)
	.def_readwrite("final_time_used",&MDT::final_time_used)
	.def_readwrite("accuracy_estimate",&MDT::accuracy_estimate)
	.def_readwrite("accuracy_estimate_user_coords",&MDT::accuracy_estimate_user_coords)
	.def_readwrite("cycle_num",&MDT::cycle_num)
	.def_readwrite("endgame_success",&MDT::endgame_success)

	.def_readwrite("function_residual",&MDT::function_residual)

	.def_readwrite("multiplicity",&MDT::multiplicity)
	.def_readwrite("is_real",&MDT::is_real)
	.def_readwrite("is_finite",&MDT::is_finite)
	.def_readwrite("is_singular",&MDT::is_singular)
	;
}




template<typename NumT>
void ExposeEndgameBoundaryMetaData(std::string const& class_name){
	using namespace bertini::algorithm;

	using MDT = EGBoundaryMetaData<NumT>;

	class_<MDT>(class_name.c_str(),init<>())
	.def_readwrite("path_point",&MDT::path_point)
	.def_readwrite("success_code",&MDT::success_code)
	.def_readwrite("last_used_stepsize",&MDT::last_used_stepsize)
	;
}


void ExposeZDMetaData(){
	using namespace bertini::algorithm;

	class_<AlgorithmMetaData>("AlgorithmMetaData",init<>())
	.def_readwrite("number_path_failures",&AlgorithmMetaData::number_path_failures)
	.def_readwrite("number_path_successes",&AlgorithmMetaData::number_path_successes)
	.def_readwrite("number_paths_tracked",&AlgorithmMetaData::number_paths_tracked)

	.def_readwrite("start_time",&AlgorithmMetaData::start_time)
	.def_readwrite("elapsed_time",&AlgorithmMetaData::elapsed_time)
	;

	ExposeSolutionMetaData<mpfr_complex>("SolutionMetaDataMultiPrec");
	ExposeSolutionMetaData<dbl_complex>("SolutionMetaDataDoublePrec");


	ExposeEndgameBoundaryMetaData<mpfr_complex>("EndgameBoundaryMetaDataMultiPrec");
	ExposeEndgameBoundaryMetaData<dbl_complex>("EndgameBoundaryMetaDataDoublePrec");

}


}} // namespaces