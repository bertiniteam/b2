// python/endgame_mp_export.cpp — MultiplePrecisionTracker endgame registrations.

#include "endgame_export.hpp"

namespace bertini{
	namespace python{

		void ExportFMPSEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedMultiplePowerSeriesEndgame",
				"The fixed but arbitrary precision implementation of the power series endgame",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings"))
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}

		void ExportFMCauchyEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FixedMultipleCauchyEndgame",
				"The fixed multiple precision implementation of the Cauchy endgame",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings"))
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

}} // namespaces
