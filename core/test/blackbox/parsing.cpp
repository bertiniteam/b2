//This file is part of Bertini 2.
//
//test/blackbox/parsing.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/blackbox/parsing.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/blackbox/parsing.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame

#include <boost/test/unit_test.hpp>

#include "bertini2/system/precon.hpp"
#include "bertini2/blackbox/global_configs.hpp"

#include "bertini2/io/parsing.hpp"


BOOST_AUTO_TEST_SUITE(blackbox_test)

BOOST_AUTO_TEST_SUITE(parsing_configs)

using namespace bertini;

BOOST_AUTO_TEST_CASE(parse1)
{
	
	using AllConfsD = blackbox::config::Configs::All<bertini::dbl>::type;
	using AllConfsMP = blackbox::config::Configs::All<bertini::mpfr>::type;

std::string config = 
R"(outputlevel: 0;
randomseed: 72;
tracktype: 1;
odepredictor: 8;
finaltol: 1e-12;
endgamenum: 2;
numsamplepoints: 8;
samplefactor: 0.8;
maxcyclenum: 6;
ampmaxprec: 3096;
ampsafetydigits1: 0;
ampsafetydigits2: 0;
maxstepsize: 0.05;
maxnumbersteps: 20000;
securitylevel: 1;
EndpointFiniteThreshold: 1e10;
pathtruncationthreshold: 1e10;
endgamebdry: 0.001;
stepsforincrease: 10;
stepfailfactor: 0.45;
stepsuccessfactor: 1.5;
functiontolerance: 1e-7;
tracktolbeforeeg: 1e-8;
tracktolduringeg: 1e-8;
sharpendigits: 60;
condnumthreshold: 1e300;
maxstepsbeforenewton: 0;
maxnewtonits: 1;)";

	auto results_double = bertini::parsing::classic::ConfigParser<AllConfsD>::Parse(config);
	auto results_mp = bertini::parsing::classic::ConfigParser<AllConfsMP>::Parse(config);
}



BOOST_AUTO_TEST_SUITE_END() // end the parsing sub-suite

BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite
