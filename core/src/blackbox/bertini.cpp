#include "bertini2/bertini.hpp"

int main(int argument_count, char** arguments)
{	
	using namespace bertini;

	auto parsed = ParseArgcArgv(argument_count, arguments);

	serial::Initialize();
	parallel::Initialize();

	int result = MainModeSwitch(parsed);

	parallel::Finalize();
	serial::Finalize();

	return result;
}
