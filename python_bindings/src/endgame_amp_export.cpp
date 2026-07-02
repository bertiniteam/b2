// python/endgame_amp_export.cpp — AMPTracker endgame registrations.

#include "endgame_export.hpp"

namespace bertini{
	namespace python{

		void ExportAMPPSEG()
		{
			using TrackerT = AMPTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("AMPPowerSeriesEndgame",
				"The adaptive precision implementation of the power series endgame.",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings") )
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}

		void ExportAMPCauchyEG()
		{
			using TrackerT = AMPTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("AMPCauchyEndgame",
				"The adaptive precision implementation of the Cauchy endgame",
				init<TrackerT const&>( (arg("self"),arg("tracker")), "Default construct with default settings"))
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

}} // namespaces
