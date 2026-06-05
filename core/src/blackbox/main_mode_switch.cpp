//This file is part of Bertini 2.
//
//src/blackbox/main_mode_switch.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/blackbox/main_mode_switch.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/blackbox/main_mode_switch.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file src/blackbox/main_mode_switch.cpp 

\brief Provides the main mode switch for the Bertini2 executable program.
*/


#include "bertini2/blackbox/main_mode_switch.hpp"
#include "bertini2/blackbox/algorithm_builder.hpp"
#include "bertini2/io/parsing/classic_utilities.hpp"
#include "bertini2/io/parsing/settings_parsers.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"

#include <fstream>
#include <iostream>

namespace bertini{

namespace {

int RunZeroDim(std::string const& config_str, std::string const& input_str)
{
	blackbox::AlgoBuilder builder;
	if (builder.ClassicBuild(config_str, input_str) != 0)
	{
		std::cerr << "error: failed to build zero-dim algorithm from input\n";
		return 1;
	}

	auto* alg = dynamic_cast<algorithm::AnyZeroDim*>(builder.GetAlg());
	if (!alg)
	{
		std::cerr << "error: algorithm builder returned wrong type\n";
		return 1;
	}

	alg->Run();

	// TODO(MPI): file writes belong on rank 0 only
	{
		std::ofstream main_data{"main_data"};
		alg->WriteMainData(main_data);
	}
	{
		std::ofstream raw_data{"raw_data"};
		alg->WriteRawData(raw_data);
	}

	std::cout << "bertini: wrote main_data and raw_data\n";
	return 0;
}

} // anonymous namespace


int MainModeSwitch(ParsedArgs const& args)
{
	// TODO(MPI): file reading belongs on rank 0 only; broadcast config_str and
	// input_str to all ranks before ClassicBuild so every rank builds an
	// identical algorithm object.
	std::string config_str, input_str;
	try {
		std::tie(config_str, input_str) = parsing::classic::SplitIntoConfigAndInput(args.input_file);
	} catch (std::exception const& e) {
		std::cerr << "error reading input file '" << args.input_file.string() << "': " << e.what() << "\n";
		return 1;
	}

	using AlgoChoice = algorithm::classic::AlgoChoice;
	auto choice = parsing::classic::FillConfigStruct<AlgoChoice>(config_str);

	switch (choice)
	{
		case AlgoChoice::ZeroDim:
			return RunZeroDim(config_str, input_str);

		case AlgoChoice::NID:
			std::cerr << "error: tracktype 1 (NID) is not yet implemented in Bertini2\n";
			return 1;

		default:
			std::cerr << "error: tracktype " << static_cast<int>(choice)
			          << " is not yet implemented in Bertini2\n";
			return 1;
	}
}

} // namespace bertini
