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
void EvalSystem(const bertini::System& S, const bertini::Vec<CType>& v, int num_times)
{
    auto J = S.Jacobian(v);
    for(int ii = 0; ii < num_times; ++ii)
    {
        S.Reset();
        J = S.Jacobian(v);
    }
}

#endif /* performance_tests_h */
