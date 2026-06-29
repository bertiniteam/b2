// python/zero_dim_double_export.cpp — DoublePrecisionTracker ZeroDim registrations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDDouble(){
			using TrackerT = bertini::tracking::DoublePrecisionTracker;
			// the start-system selector enum + factory are not tracker-specific; register them once,
			// here, and BEFORE the ZeroDim classes (whose 2nd constructor takes a StartSystemFactory).
			ExportStartSystemEnum();
			ExportStartSystemFactory();
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("ZeroDimSolverPowerSeriesDoublePrecision");
			ExportZeroDimSpecific<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("ZeroDimSolverCauchyDoublePrecision");

			// the User start system is not tracker-specific; register it once, here.
			ExportUserStartSystem();
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG>("HomotopySolverPowerSeriesDoublePrecision");
			ExportZeroDimUserHomotopy<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy>("HomotopySolverCauchyDoublePrecision");
		}

}} // namespaces
