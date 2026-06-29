// python/zero_dim_mp_export.cpp — MultiplePrecisionTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDMP(){
			using TrackerT = bertini::tracking::MultiplePrecisionTracker;
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimPowerSeriesFixedMultiplePrecision");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimCauchyFixedMultiplePrecision");

			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimPowerSeriesFixedMultiplePrecisionUserHomotopy");
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimCauchyFixedMultiplePrecisionUserHomotopy");
		}

}} // namespaces
