
#pragma once

#include <bertini2/nag_algorithms/zero_dim_solve.hpp>
#include <bertini2/nag_algorithms/output.hpp>
#include <bertini2/endgames.hpp>
#include <bertini2/system.hpp>

namespace demo{


using TrackerT = bertini::tracking::AMPTracker;
using Tolerances = bertini::algorithm::TolerancesConfig;
using EndgameConfT = bertini::endgame::EndgameConfig;

auto StepOne(bertini::System const& sys)
{
	using namespace bertini;
	using namespace algorithm;

	using EndgameT = typename endgame::EndgameSelector<TrackerT>::Cauchy;
	

	auto zd = bertini::algorithm::ZeroDim<TrackerT, EndgameT, bertini::System, bertini::start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-5;
	tols.newton_during_endgame = 1e-6;
	zd.Set(tols);

	auto& tr = zd.GetTracker();

	tr.SetPredictor(bertini::tracking::Predictor::HeunEuler);
	tracking::GoryDetailLogger<TrackerT> tr_logger;
	// tr.AddObserver(&tr_logger);

	endgame::GoryDetailLogger<EndgameT> eg_logger;
	zd.GetEndgame().AddObserver(&eg_logger);

	auto eg = zd.GetFromEndgame<EndgameConfT>();
	eg.final_tolerance = 1e-11;
	zd.SetToEndgame(eg);

	zd.Solve();

	return output::NonsingularSolutions::Extract(zd);
}
	
template <typename SolnContT>
auto StepTwo(bertini::System const& target_sys, bertini::System const& start_sys, bertini::System const& homotopy, SolnContT const& solns)
{
	using namespace bertini;
	using namespace tracking;
	using namespace algorithm;

	auto userss = bertini::start_system::User(start_sys, solns);

	auto zd = bertini::algorithm::ZeroDim<TrackerT, typename bertini::endgame::EndgameSelector<TrackerT>::Cauchy, bertini::System, bertini::start_system::User, bertini::policy::RefToGiven>(target_sys, userss, homotopy);

	zd.DefaultSetup();
	
	zd.GetTracker().SetPredictor(bertini::tracking::Predictor::HeunEuler);

	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-6;
	tols.newton_during_endgame = 1e-7;
	zd.Set(tols);

	auto eg = zd.GetFromEndgame<EndgameConfT>();
	eg.final_tolerance = 1e-12;
	zd.SetToEndgame(eg);
	
	zd.Solve();

	return output::AllSolutions::Extract(zd);
}


} // namespace demo

