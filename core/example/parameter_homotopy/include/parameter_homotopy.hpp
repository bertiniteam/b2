
#pragma once

#include <bertini2/nag_algorithms/zero_dim_solve.hpp>
#include <bertini2/endgames.hpp>
#include <bertini2/system.hpp>

namespace demo{


using TrackerT = bertini::tracking::AMPTracker;
using Tolerances = bertini::algorithm::TolerancesConfig;
using EndgameConfT = bertini::endgame::EndgameConfig;

// not really supposed to be void
void StepOne(bertini::System const& sys)
{
	using namespace bertini;
	using namespace tracking;
	using namespace algorithm;

	
	

	auto zd = bertini::algorithm::ZeroDim<TrackerT, typename endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-5;
	tols.newton_during_endgame = 1e-6;
	zd.Set(tols);

	auto& tr = zd.GetTracker();

	tr.SetPredictor(bertini::tracking::Predictor::HeunEuler);
	GoryDetailLogger<TrackerT> logger;
	// tr.AddObserver(&logger);


	auto eg = zd.GetFromEndgame<EndgameConfT>();
	eg.final_tolerance = 1e-11;
	zd.SetToEndgame(eg);

	zd.Solve();

	std::cout << "solved start system\n";
}

void StepTwo(bertini::System const& homotopy, bertini::System const& start_sys)
{


}


} // namespace demo

