#include "parameter_homotopy.hpp"
#include "random.hpp"
#include "my_system.hpp"

#include <chrono>
int main()
{
    bertini::LoggingInit();


    

    auto step1_params = demo::MakeStep1Parameters();
    auto target_sys_step1 = demo::ConstructSystem(step1_params);

    std::cout << "your target_system for step 1:\n\n";
    std::cout << target_sys_step1 << '\n';

    std::cout << "\n\nwith parameter values:\n\n";
    for (const auto& p : step1_params)
        std::cout << p << " ";
    std::cout << '\n';

    // now to solve the start system.
    auto stepone_solutions = demo::StepOne(target_sys_step1);

    std::cout << "done computing the " << stepone_solutions.size() << " step1 solutions, and here they are: \n";
    for (auto& iter : stepone_solutions)
        std::cout << iter << '\n' << '\n';


    auto t = bertini::MakeVariable("t");
    auto step2_stuff = demo::MakeStep2Parameters(step1_params, t);
    auto homotopy_sys_step2 = demo::ConstructSystem(std::get<0>(step2_stuff));
    homotopy_sys_step2.AddPathVariable(t);
    auto target_sys_step2 = demo::ConstructSystem(std::get<1>(step2_stuff));
    
    int num_to_step2s = 100;
    auto start_allstep2 = std::chrono::high_resolution_clock::now();
    for (int ii=0; ii<num_to_step2s; ii++)
    {
        auto start_iteration = std::chrono::high_resolution_clock::now();

        bertini::DefaultPrecision(30);
        // iterate over the parameter values.  set the 
        for (auto& p : std::get<1>(step2_stuff))
        {
            bertini::mpfr v;
            bertini::RandomReal(v, 30);
            p->precision(30);
            p->set_current_value(bertini::dbl(v));
            p->set_current_value(v);
        }

        bertini::DefaultPrecision(16);

        // std::cout << "your target system for step2:\n" << target_sys_step2 << '\n';
        // for (auto& p : std::get<1>(step2_stuff))
        //     std::cout << "solving for parameter values " << *p <<  " " << p->Eval<bertini::dbl>() << '\n';

        auto steptwo_solutions = demo::StepTwo(target_sys_step2, target_sys_step1, homotopy_sys_step2, stepone_solutions);
        
        // std::cout << "done computing the " << steptwo_solutions.size() << " step2 solutions, and here they are: \n";
        // for (auto& iter : steptwo_solutions)
        //     std::cout << iter << '\n' << '\n';

        std::cout << (std::chrono::high_resolution_clock::now() - start_iteration).count() << '\n';
    }
    std::cout << "solving " << num_to_step2s << " " << (std::chrono::high_resolution_clock::now() - start_allstep2).count() << '\n';

	return 0;
}