//
//  performance_tests.hpp
//  Xcode_b2
//
//  Created by Jeb Collins University of Mary Washington. All rights reserved.
//

#ifndef performance_tests_h
#define performance_tests_h
#include "construct_system.hpp"

template<typename CType>
void EvalAndLUTests(const bertini::System& S, const bertini::Vec<CType>& v,
                      const bertini::Vec<CType>& b, int num_times)
{
    auto J = S.Jacobian(v);
    for(int ii = 0; ii < num_times; ++ii)
    {
        S.Reset();
        J = S.Jacobian(v);
        J.lu().solve(b);
    }
}

template<typename CType>
void MatrixMultTests(const Mat<CType>& A, int num_times)
{
    for(int ii = 0; ii < num_times; ++ii)
    {
        A*A;
    }
}

#endif /* performance_tests_h */
