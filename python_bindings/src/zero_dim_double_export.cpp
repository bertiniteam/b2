// python/zero_dim_double_export.cpp — DoublePrecisionTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDDouble(){
			using TrackerT = bertini::tracking::DoublePrecisionTracker;
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesDoublePrecisionTotalDegree");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyDoublePrecisionTotalDegree");
		}

}} // namespaces
