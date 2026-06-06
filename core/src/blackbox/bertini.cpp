#include "bertini2/bertini.hpp"

int main(int argument_count, char** arguments)
{	
	using namespace bertini;

	auto parsed = ParseArgcArgv(argument_count, arguments);

	parallel::Initialize();  // MPI_Init first so rank is known
	serial::Initialize();    // splash on rank 0 only

	int result = MainModeSwitch(parsed);

	parallel::Finalize();
	serial::Finalize();

	return result;
}
