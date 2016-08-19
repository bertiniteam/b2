#include "bertini2/bertini.hpp"

int main(int argument_count, char** arguments)
{	
	using namespace bertini;

	ParseArgcArgv(argument_count, arguments);

	serial::Initialize();
	parallel::Initialize();



	MainModeSwitch();



	parallel::Finalize();
	serial::Finalize();

	return 0;
}
