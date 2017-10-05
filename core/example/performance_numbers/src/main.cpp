#include "performance_tests.hpp"
#include <ctime>


int main()
{

    int num_evaluations = 1; ///> Number of times to evaluate the Jacobian in each run
    int num_test_runs = 5000; ///> number of times to run the test for average
    int num_precisions = 10;  ///> number of different precisions to use
    int max_precision = 308; ///> maximum precision used for testing
    int matrix_N = 100; ///> size of matrix for matrix multiplication
    
    std::vector<int> precisions(num_precisions-1);
    for(int P = 0; P < num_precisions-1; ++P)
    {
        precisions[P] = std::floor(16 + ((max_precision)-16.0)/num_precisions*(P));
    }
    precisions.push_back(max_precision);
    
    
    // Compute the time using CPU clock time, not wall clock time.
    auto start = std::clock();
    auto end = std::clock();

    auto sys1 = demo::ConstructSystem1();

    
    
    
    
    auto v_d = demo::GenerateSystemInput<dbl>(sys1);
    auto b_d = demo::GenerateRHS<dbl>(sys1);
    auto A_d = demo::GenerateMatrix<dbl>(matrix_N);
    auto v_mp = demo::GenerateSystemInput<mpfr>(sys1);
    auto b_mp = demo::GenerateRHS<mpfr>(sys1);
    auto A_mp = demo::GenerateMatrix<mpfr>(matrix_N);
    
    
    //Get base number for double precision
    std::cout << "\n\n\nTesting system Jacobian evaluation...:\n\n";
    double time_delta_d = 0;
    for(int ii = 0; ii < num_test_runs; ++ii)
    {
        start = std::clock();
        EvalAndLUTests(sys1, v_d, b_d, num_evaluations);
        MatrixMultTests(A_d, 100*num_evaluations);
        end = std::clock();
        time_delta_d += (double)(end-start)/(double)(CLOCKS_PER_SEC);
    }
    time_delta_d = time_delta_d/num_test_runs;

    
    std::cout << "Average time taken:\n\n";
    std::cout << time_delta_d << std::endl;

    
    // Now work with various precisions for mpfr
    std::cout << "Evaluating Jacobian in multiple precision:\n\n";
    Vec<double> time_delta_mp(num_precisions);
    for(int PP = 0; PP < num_precisions; ++PP)
    {
        std::cout << "Evaluating with precision " << precisions[PP] << "...\n";
        time_delta_mp(PP) = 0;
        v_mp = demo::GenerateSystemInput<mpfr>(sys1, precisions[PP]);
        b_mp = demo::GenerateRHS<mpfr>(sys1, Precision(v_mp));
        sys1.precision(Precision(v_mp));
        for(int ii = 0; ii < num_test_runs; ++ii)
        {
            start = std::clock();
            EvalAndLUTests(sys1, v_mp, b_mp, num_evaluations);
            MatrixMultTests(A_mp, 100*num_evaluations);
            end = std::clock();
            time_delta_mp(PP) += (double)(end-start)/(double)(CLOCKS_PER_SEC);
        }
        time_delta_mp(PP) = time_delta_mp(PP)/num_test_runs;
    }
    
    
    
//    std::cout << time_delta_mp << std::endl;
    
    auto time_factors = time_delta_mp/time_delta_d;
    
    
    // Compute coefficient for linear fit
    Mat<double> M(2,2);
    Vec<double> b(2);
    
    M(0,0) = num_precisions;
    M(0,1) = 0;
    M(1,1) = 0;
    b(0) = 0; b(1) = 0;
    for(int ii = 0; ii < num_precisions; ++ii)
    {
        M(0,1) += precisions[ii];
        M(1,1) += pow(precisions[ii],2);
        b(0) += time_factors(ii);
        b(1) += precisions[ii]*time_factors(ii);
    }
    M(1,0) = M(0,1);
    
    Vec<double> x = M.lu().solve(b);
    
//    std::cout << x(0) << std::endl;
    std::cout << "y(P) = "<< x(1)<<"x + "<< x(0) << std::endl;
    
    
	return 0;
}
