#include "parameter_homotopy.hpp"
#include "random.hpp"
#include "my_system.hpp"


int main()
{
    bertini::LoggingInit();



    auto step1_params = demo::MakeStep1Parameters();

    auto target_sys = demo::ConstructSystem(step1_params);
    auto start_sys = demo::ConstructStart(target_sys);

	auto homotopy = demo::ConstructHomotopy(target_sys, start_sys);


	std::cout << "your homotopy:\n\n";
	std::cout << homotopy << '\n';

    std::cout << "your start_sys:\n\n";
    std::cout << start_sys << '\n';

    std::cout << "your target_sys:\n\n";
    std::cout << target_sys << '\n';


    // now to solve the start system.
    auto stepone_solutions = demo::StepOne(start_sys);

    demo::StepTwo(target_sys, start_sys, homotopy, stepone_solutions);
    
	return 0;
}