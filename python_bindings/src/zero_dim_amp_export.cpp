// python/zero_dim_amp_export.cpp — AMPTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDAMP(){
			using TrackerT = bertini::tracking::AMPTracker;
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesAdaptivePrecisionTotalDegree");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyAdaptivePrecisionTotalDegree");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::MHomogeneous>("ZeroDimPowerSeriesAdaptivePrecisionMHomogeneous");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::MHomogeneous>("ZeroDimCauchyAdaptivePrecisionMHomogeneous");

			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimPowerSeriesAdaptivePrecisionUserHomotopy");
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimCauchyAdaptivePrecisionUserHomotopy");
		}

}} // namespaces
