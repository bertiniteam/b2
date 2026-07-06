// python/zero_dim_amp_export.cpp — AMPTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDAMP(){
			using TrackerT = bertini::tracking::AMPTracker;
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimSolverPowerSeriesAdaptivePrecision");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimSolverCauchyAdaptivePrecision");

			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("HomotopySolverPowerSeriesAdaptivePrecision");
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("HomotopySolverCauchyAdaptivePrecision");
		}

}} // namespaces
