// python/endgame_double_export.cpp — DoublePrecisionTracker endgame registrations.

#include "endgame_export.hpp"

namespace bertini{
	namespace python{

		void ExportFDPSEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedDoublePSEG",
				"The double-precision implementation of the power series endgame",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings"))
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}

		void ExportFDCauchyEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FixedDoubleCauchyEG",
				"The fixed double precision implementation of the Cauchy endgame",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings"))
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

}} // namespaces
