#include "performance_tests.hpp"
#include <ctime>


int main()
{


    

    auto sys1 = demo::ConstructSystem1();

    std::cout << "\n\n\nTesting system Jacobian evaluation...:\n\n";
    
    
    sys1.Differentiate();
    
    bertini::Vec<dbl> v_d(4);
    v_d << dbl(2.435, -1.6748),dbl(-1.34, 5.231850),dbl(-6.3728467, 2.89570486),dbl(3.4957219562, -4.098763882);
    bertini::Vec<mpfr> v_mp(4);
    v_mp << mpfr("2.435", "-1.6748"),mpfr("-1.34", "5.231850"),mpfr("-6.3728467", "2.89570486"),mpfr("3.4957219562", "-4.098763882");

    
    
    auto start = std::clock();
    sys1.precision(16);
    EvalSystem(sys1, v_mp, 100);
    auto end = std::clock();
    
    std::cout << "time taken:\n\n";
    std::cout << end-start << std::endl;

//    std::cout << "\n\nwith parameter values:\n\n";
//    for (const auto& p : step1_params)
//        std::cout << p << " ";
//    std::cout << '\n';
//
//    // now to solve the start system.
//    auto stepone_solutions = demo::StepOne(target_sys_step1);
//
//    std::cout << "done computing the " << stepone_solutions.size() << " step1 solutions, and here they are: \n";
//    for (auto& iter : stepone_solutions)
//        std::cout << iter << '\n' << '\n';
//
//
//    auto t = bertini::MakeVariable("t");
//    auto step2_stuff = demo::MakeStep2Parameters(step1_params, t);
//    auto homotopy_sys_step2 = demo::ConstructSystem(std::get<0>(step2_stuff));
//    homotopy_sys_step2.AddPathVariable(t);
//    auto target_sys_step2 = demo::ConstructSystem(std::get<1>(step2_stuff));
//
//
//
//    bertini::DefaultPrecision(30);
//    // iterate over the parameter values.  set the
//    for (auto& p : std::get<1>(step2_stuff))
//    {
//        bertini::mpfr v;
//        bertini::RandomReal(v, 30);
//        p->precision(30);
//        p->set_current_value(bertini::dbl(v));
//        p->set_current_value(v);
//    }
//
//    bertini::DefaultPrecision(16);
//
//    std::cout << "your target system for step2:\n" << target_sys_step2 << '\n';
//    for (auto& p : std::get<1>(step2_stuff))
//        std::cout << "solving for parameter values " << *p <<  " " << p->Eval<bertini::dbl>() << '\n';
//
//
//
//    auto steptwo_solutions = demo::StepTwo(target_sys_step2, target_sys_step1, homotopy_sys_step2, stepone_solutions);
//
//    std::cout << "done computing the " << steptwo_solutions.size() << " step2 solutions, and here they are: \n";
//    for (auto& iter : steptwo_solutions)
//        std::cout << iter << '\n' << '\n';

	return 0;
}
