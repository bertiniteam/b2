// python/nid_amp_export.cpp — AMPTracker NID algorithm registrations.

#include "numerical_irreducible_decomposition_export.hpp"

namespace bertini{
	namespace python{

		void ExportNIDAMP(){
			using TrackerT = bertini::tracking::AMPTracker;
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("NIDPowerSeriesAdaptivePrecision");
			ExportNIDSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("NIDCauchyAdaptivePrecision");
		}

}} // namespaces
