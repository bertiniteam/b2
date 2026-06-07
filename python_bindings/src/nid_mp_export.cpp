// python/nid_mp_export.cpp — MultiplePrecisionTracker NID algorithm registrations.

#include "numerical_irreducible_decomposition_export.hpp"

namespace bertini{
	namespace python{

		void ExportNIDMP(){
			using TrackerT = bertini::tracking::MultiplePrecisionTracker;
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesFixedMultiplePrecision");
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyFixedMultiplePrecision");
		}

}} // namespaces
