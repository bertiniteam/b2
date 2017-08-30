//This file is part of Bertini 2.
//
//zero_dim.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//zero_dim.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with zero_dim.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

/**
\file test/nag_algorithms/zero_dim.cpp  Tests the zero dim algorithm.
*/

// individual authors of this file include:
// dani brake, university of notre dame

#include "bertini2/system/precon.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"
#include <boost/test/unit_test.hpp>
#include "bertini2/nag_algorithms/output.hpp"


BOOST_AUTO_TEST_SUITE(zero_dim)

using TrackerT = bertini::tracking::DoublePrecisionTracker;



BOOST_AUTO_TEST_CASE(can_run_griewank_osborn)
{
	using namespace bertini;
	using namespace tracking;

	using Tolerances = algorithm::TolerancesConfig;
	using EndgameConfT = endgame::EndgameConfig;
	

	auto sys = system::Precon::GriewankOsborn();

	auto zd = algorithm::ZeroDim<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::Cauchy, decltype(sys), start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-5;
	tols.newton_during_endgame = 1e-6;
	zd.Set(tols);

	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(&logger);


	auto eg = zd.GetFromEndgame<EndgameConfT>();
	eg.final_tolerance = 1e-12;
	zd.SetToEndgame(eg);

	zd.Solve();


	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}




BOOST_AUTO_TEST_CASE(can_run_change_some_settings)
{
	
	using namespace bertini;
	using namespace tracking;

	using Tolerances = algorithm::TolerancesConfig;

	auto sys = system::Precon::GriewankOsborn();

	auto zd = algorithm::ZeroDim<TrackerT, bertini::endgame::EndgameSelector<TrackerT>::PSEG, decltype(sys), start_system::TotalDegree>(sys);

	zd.DefaultSetup();
	
	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(&logger);

	

	auto tols = zd.Get<Tolerances>();
	tols.newton_before_endgame = 1e-4;
	tols.newton_during_endgame = 1e-4;

	zd.Set(tols);

	zd.Solve();


	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}



/**
Check whether we can run zero dim on the non-homogenized version of Griewank Osborn.
*/
BOOST_AUTO_TEST_CASE(reference_managed_systems_GO_nonhom)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	auto TD = start_system::TotalDegree(sys);

	auto t = MakeVariable("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::ZeroDim<
				TrackerT, 
				bertini::endgame::EndgameSelector<TrackerT>::PSEG, 
				decltype(sys), 
				start_system::TotalDegree,
				policy::RefToGiven
					>
			(sys, TD, h); 
	// have to pass in all three to the constructor, because using references. 

	zd.DefaultSetup();

	zd.Solve();

	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}





/**
Check whether we can run zero dim on the non-homogenized version of Griewank Osborn.
*/
BOOST_AUTO_TEST_CASE(reference_managed_systems_GO)
{
	using namespace bertini;
	using namespace tracking;

	auto sys = system::Precon::GriewankOsborn();
	sys.Homogenize();
	sys.AutoPatch();

	auto TD = start_system::TotalDegree(sys);

	auto t = MakeVariable("t");
	auto h = (1-t)* sys + t*TD;
	h.AddPathVariable(t);

	auto zd = algorithm::ZeroDim<
				TrackerT, 
				bertini::endgame::EndgameSelector<TrackerT>::PSEG, 
				decltype(sys), 
				start_system::TotalDegree,
				policy::RefToGiven
					>
			(sys, TD, h); 
	// have to pass in all three to the constructor, because using references. 

	zd.DefaultSetup();

	auto& tr = zd.GetTracker();
	GoryDetailLogger<TrackerT> logger;
	tr.AddObserver(&logger);
	
	zd.Solve();

	bertini::algorithm::output::Classic<decltype(zd)>::All(std::cout, zd);
}


BOOST_AUTO_TEST_SUITE_END()
