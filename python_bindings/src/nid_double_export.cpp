// python/nid_double_export.cpp — DoublePrecisionTracker NID algorithm registrations.

#include "numerical_irreducible_decomposition_export.hpp"

namespace bertini{
	namespace python{

		void ExportNIDDouble(){
			using TrackerT = bertini::tracking::DoublePrecisionTracker;
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesDoublePrecision");
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyDoublePrecision");
		}

}} // namespaces
