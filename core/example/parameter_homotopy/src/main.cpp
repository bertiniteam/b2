#include "parameter_homotopy.hpp"
#include "random.hpp"
#include "my_system.hpp"


int main()
{
    bertini::LoggingInit();


    auto t = bertini::MakeVariable("t");

    auto step1_params = demo::MakeStep1Parameters();

    auto step2_stuff = demo::MakeStep2Parameters(step1_params, t);


    auto target_sys_step1 = demo::ConstructSystem(step1_params);

    std::cout << "your target_system for step 1:\n\n";
    std::cout << target_sys_step1 << '\n';

    // now to solve the start system.
    auto stepone_solutions = demo::StepOne(target_sys_step1);

    std::cout << "done computing step1 solutions, and here they are: \n";
    for (auto& iter : stepone_solutions)
        std::cout << iter << '\n';


    auto homotopy_sys_step2 = demo::ConstructSystem(std::get<0>(step2_stuff));
    homotopy_sys_step2.AddPathVariable(t);
    auto target_sys_step2 = demo::ConstructSystem(std::get<1>(step2_stuff));
    


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

    std::cout << "your target system for step2:\n" << target_sys_step2 << '\n';



    for (auto& p : std::get<1>(step2_stuff))
        std::cout << "solving for parameter values " << *p <<  " " << p->Eval<bertini::dbl>() << '\n';

    auto steptwo_solutions = demo::StepTwo(target_sys_step2, target_sys_step1, homotopy_sys_step2, stepone_solutions);
    
    std::cout << "done computing step2 solutions, and here they are: \n";
    for (auto& iter : steptwo_solutions)
        std::cout << iter << '\n';

	return 0;
}