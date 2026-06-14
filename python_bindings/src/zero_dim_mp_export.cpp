// python/zero_dim_mp_export.cpp — MultiplePrecisionTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDMP(){
			using TrackerT = bertini::tracking::MultiplePrecisionTracker;
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::TotalDegree>("ZeroDimPowerSeriesFixedMultiplePrecisionTotalDegree");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>("ZeroDimCauchyFixedMultiplePrecisionTotalDegree");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, bertini::System, bertini::start_system::MHomogeneous>("ZeroDimPowerSeriesFixedMultiplePrecisionMHomogeneous");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::MHomogeneous>("ZeroDimCauchyFixedMultiplePrecisionMHomogeneous");

			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimPowerSeriesFixedMultiplePrecisionUserHomotopy");
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimCauchyFixedMultiplePrecisionUserHomotopy");
		}

}} // namespaces
