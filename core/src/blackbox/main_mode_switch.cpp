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
#include "bertini2/parallel.hpp"
#include "bertini2/random.hpp"

#include <fstream>
#include <iostream>

#ifdef BERTINI2_HAVE_MPI
#include "bertini2/parallel/mpi_include.hpp"
#endif


namespace bertini{

namespace {

int RunZeroDim(std::string const& config_str, std::string const& input_str)
{
	// Parse and apply the RNG seed before any setup draws (gamma, TD-constants, patch).
	auto rand_cfg = parsing::classic::FillConfigStruct<algorithm::RandomConfig>(config_str);

#ifdef BERTINI2_HAVE_MPI
	// Manager sets the seed (possibly from entropy), then broadcasts the effective
	// (non-zero) seed to workers so all ranks share identical per-path streams.
	if (parallel::IsManager())
		SetGlobalSeed(rand_cfg.random_seed);
	unsigned long effective_seed = GetGlobalSeed();
	MPI_Bcast(&effective_seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if (!parallel::IsManager())
		SetGlobalSeed(effective_seed);
#else
	SetGlobalSeed(rand_cfg.random_seed);
#endif

	std::cout << "bertini: random seed = " << GetGlobalSeed() << "\n";

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

	// The structured output directory is simply part of the program's output (like
	// main_data): produced always, freely deletable, no flag.  BERTINI_RECORDS_DIR
	// overrides the default location; setting it EMPTY is the off switch.  Manager
	// rank only; workers never touch records.
	if (parallel::IsManager())
	{
		char const* records_dir = std::getenv("BERTINI_RECORDS_DIR");
		if (records_dir && *records_dir == '\0')
			std::cout << "bertini: records off (BERTINI_RECORDS_DIR is empty)\n";
		else
		{
			alg->RecordToPath(records_dir ? records_dir : "bertini_output");
			std::cout << "bertini: records at " << (records_dir ? records_dir : "bertini_output") << "\n";
		}
	}

	alg->Run();

	if (parallel::IsManager())
	{
		// archive the ACTUAL input file in the records: the truest statement of what was
		// asked, referenced from this run as a `given`
		{
			auto const input_definition =
				alg->PutRecordsDefinition(config_str + input_str, "givens", "cli_input");
			if (!input_definition.empty())
				alg->AppendRecordJson(std::string("{\"kind\":\"given\",\"role\":\"cli_input\",\"run\":\"")
					+ alg->RecordsRunIdentity() + "\",\"source\":\"" + input_definition + "\"}");
		}
		{
			std::ofstream main_data{"main_data"};
			alg->WriteMainData(main_data);
		}
		{
			std::ofstream raw_data{"raw_data"};
			alg->WriteRawData(raw_data);
		}
		// Bertini 1.7-compatible machine-readable solution files (count-led coordinate blocks),
		// so tooling that parses Bertini 1 output reads these unchanged.
		{
			std::ofstream f{"finite_solutions"};
			alg->WriteFiniteSolutions(f);
		}
		{
			std::ofstream f{"real_finite_solutions"};
			alg->WriteRealFiniteSolutions(f);
		}
		{
			std::ofstream f{"nonsingular_solutions"};
			alg->WriteNonsingularSolutions(f);
		}
		{
			std::ofstream f{"singular_solutions"};
			alg->WriteSingularSolutions(f);
		}
		{
			std::ofstream f{"raw_solutions"};
			alg->WriteRawSolutions(f);
		}
		std::cout << "bertini: wrote main_data, raw_data, and solution files\n";
	}
	return 0;
}

} // anonymous namespace


int MainModeSwitch(ParsedArgs const& args)
{
	std::string config_str, input_str;

#ifdef BERTINI2_HAVE_MPI
	// Rank 0 reads the input file, then broadcasts both strings to all ranks so
	// every rank builds an identical algorithm object from the same source text.
	if (parallel::IsManager())
	{
		try {
			std::tie(config_str, input_str) = parsing::classic::SplitIntoConfigAndInput(args.input_file);
		} catch (std::exception const& e) {
			std::cerr << "error reading input file '" << args.input_file.string() << "': " << e.what() << "\n";
			// Broadcast empty strings so workers don't hang waiting for broadcast.
			parallel::mpi_broadcast_string(parallel::WorldComm(), config_str, 0);
			parallel::mpi_broadcast_string(parallel::WorldComm(), input_str,  0);
			return 1;
		}
	}
	parallel::mpi_broadcast_string(parallel::WorldComm(), config_str, 0);
	parallel::mpi_broadcast_string(parallel::WorldComm(), input_str,  0);
	if (config_str.empty() && input_str.empty())
		return 1;  // rank 0 failed to read; workers bail out cleanly
#else
	try {
		std::tie(config_str, input_str) = parsing::classic::SplitIntoConfigAndInput(args.input_file);
	} catch (std::exception const& e) {
		std::cerr << "error reading input file '" << args.input_file.string() << "': " << e.what() << "\n";
		return 1;
	}
#endif

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
