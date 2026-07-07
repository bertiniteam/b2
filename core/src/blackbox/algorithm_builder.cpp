//This file is part of Bertini 2.
//
//src/blackbox/algorithm_builder.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/blackbox/algorithm_builder.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/blackbox/algorithm_builder.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file src/blackbox/algorithm_builder.cpp 

\brief Provides the methods for building algorithms from files or streamable sources.
*/

#include "bertini2/blackbox/algorithm_builder.hpp"
#include "bertini2/blackbox/global_configs.hpp"
#include "bertini2/blackbox/switches_zerodim.hpp"
#include "bertini2/io/parsing/settings_parsers.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"

#include <iostream>

namespace bertini{
namespace blackbox{

int AlgoBuilder::ClassicBuild(std::string const& config_str, std::string const& input_str)
{
	// Parse the polynomial system from the INPUT section
	System sys;
	try {
		sys = System{input_str};
	} catch (std::exception const& e) {
		std::cerr << "error: failed to parse system from input section: " << e.what() << "\n";
		return 1;
	}

	// Parse all configuration structs (double precision versions suffice for
	// choosing algorithm types; mpfr versions would be used for mp-specific defaults)
	using AllConfsD = config::Configs::All<complex_dbl>::type;
	decltype(parsing::classic::ConfigParser<AllConfsD>::Parse(config_str)) cfgs_d;
	try {
		cfgs_d = parsing::classic::ConfigParser<AllConfsD>::Parse(config_str);
	} catch (std::exception const& e) {
		std::cerr << "error: failed to parse configuration: " << e.what() << "\n";
		return 1;
	}

	// Select tracker type from PrecisionType (mptype in Bertini1 syntax):
	//   Fixed         (mptype: 0) -> FixedDouble
	//   FixedMultiple (mptype: 1) -> FixedMultiple (fixed-precision multi)
	//   Adaptive      (mptype: 2, default) -> Adaptive (AMP)
	auto prec_type = std::get<tracking::PrecisionType>(cfgs_d);
	type::Tracker tracker_type;
	switch (prec_type) {
		case tracking::PrecisionType::Fixed:          tracker_type = type::Tracker::FixedDouble;   break;
		case tracking::PrecisionType::FixedMultiple:  tracker_type = type::Tracker::FixedMultiple; break;
		default:                                      tracker_type = type::Tracker::Adaptive;       break;
	}

	// Select endgame type from endgamenum (Bertini1 syntax):
	//   1 -> PowerSeries (PSEG), 2 -> Cauchy (default)
	using algorithm::classic::EndgameChoiceConfig;
	using algorithm::classic::EndgameChoice;
	auto eg_choice = std::get<EndgameChoiceConfig>(cfgs_d).endgame;
	type::Endgame endgame_type = (eg_choice == EndgameChoice::PowerSeries)
	                              ? type::Endgame::PowerSeries
	                              : type::Endgame::Cauchy;

	// Infer the start system from the variable-group structure, the way classic
	// Bertini does: a single affine variable group -> total degree; multiple
	// variable groups or any homogeneous variable group -> multihomogeneous.
	// (User-defined homotopies still need a dedicated input section.)
	type::Start start_type = InferStartType(sys);

	ZeroDimRT rt{start_type, tracker_type, endgame_type};

	std::unique_ptr<algorithm::AnyZeroDim> zd_alg;
	try {
		zd_alg = MakeZeroDim(rt, sys);
	} catch (std::exception const& e) {
		std::cerr << "error: failed to instantiate zero-dim algorithm: " << e.what() << "\n";
		return 1;
	}

	zd_alg->ApplyParsedConfigs(config_str);
	alg_ = std::move(zd_alg);

	return 0;
}

} // namespace blackbox
} // namespace bertini

