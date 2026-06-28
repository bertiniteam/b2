#include "bertini2/bertini.hpp"
#include "bertini2/fast_allocator.hpp"

int main(int argument_count, char** arguments)
{
	using namespace bertini;

	// Route GMP/MPFR/MPC limb allocation through mimalloc (if built with BERTINI2_FAST_ALLOC).
	// First thing, before any multiprecision work.  No-op if disabled.
	InstallFastAllocator();

	auto parsed = ParseArgcArgv(argument_count, arguments);

	parallel::Initialize();  // MPI_Init first so rank is known
	serial::Initialize();    // splash on rank 0 only

	int result = MainModeSwitch(parsed);

	parallel::Finalize();
	serial::Finalize();

	return result;
}
