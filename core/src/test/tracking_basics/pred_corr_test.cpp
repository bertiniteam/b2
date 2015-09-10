//This file is part of Bertini 2.0.
//
//start_system_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  start_system_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

//start_system_test.cpp
//


#include <boost/test/unit_test.hpp>
#include "start_system.hpp"

#include <boost/multiprecision/mpfr.hpp>

#include "mpfr_complex.hpp"
#include "tracking.hpp"


using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = boost::multiprecision::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(tracking_basics)


BOOST_AUTO_TEST_CASE(circle_line_euler_1_corr_step)
{
    boost::multiprecision::mpfr_float::default_precision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);


    // Starting point in spacetime step
    dbl xn_d(2.3,0.2);
    dbl yn_d(1.1, 1.87);
    mpfr xn_mp("2.3","0.2");
    mpfr yn_mp("1.1", "1.87");
    // Starting time
    dbl current_time = .9;
    // Time step
    dbl delta_t = .1;
    
    
    
    bertini::System sys;
    Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
    
    VariableGroup vars{x,y};
    
    sys.AddVariableGroup(vars);
    sys.AddPathVariable(t);
    
    // Define homotopy system
    sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
    sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
    
    
    
    auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

    
    mpfr_float tracking_tolerance_mp("1e-5");
    double tracking_tolerance_d(tracking_tolerance_mp);


    // Exact value of step after Euler prediction
    Vec<dbl> predicted_d(2);
    predicted_d << dbl(2.19163044775416716032722848776575,-0.203532362690064835307109968099797),
        dbl(1.82976707889226183423896735889061, 2.37621890895572354790101570972854);

    Vec<mpfr> predicted_mp(2);
    predicted_mp << mpfr("2.19163044775416716032722848776575","-0.203532362690064835307109968099797"),
        mpfr("1.82976707889226183423896735889061", "2.37621890895572354790101570972854");
    

    // Exact value of step after a one-Newton-step correction
    Vec<dbl> corrected_d(2);
    corrected_d << dbl(1.27914921783237417760452216130584,0.287985101426454933172731469315295),
        dbl(0.16018906270391684942121729748759, -0.0639966892058788740383847709589545);

    Vec<mpfr> corrected_mp(2);
    corrected_mp << mpfr("1.27914921783237417760452216130584","0.287985101426454933172731469315295"),
        mpfr("0.16018906270391684942121729748759", "-0.0639966892058788740383847709589545");

    


    
    
}



BOOST_AUTO_TEST_CASE(circle_line_euler_2_corr_step)
{
    boost::multiprecision::mpfr_float::default_precision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);

    // Starting point in spacetime step
    dbl xn_d(2.3,0.2);
    dbl yn_d(1.1, 1.87);
    mpfr xn_mp("2.3","0.2");
    mpfr yn_mp("1.1", "1.87");
    // Starting time

    Vec<dbl> current_space(2);
    current_space << xn_d, yn_d;

    dbl current_time = .9;
    // Time step
    dbl delta_t = .1;

    
    
    
    bertini::System sys;
    Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
    
    VariableGroup vars{x,y};
    
    sys.AddVariableGroup(vars);
    sys.AddPathVariable(t);
    
    // Define homotopy system
    // [t(x^2-1) + (1-t)(x^2 + y^2 -4);...
    //   t(y-1) + (1-t)(2x+5y)]
    sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
    sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
    
    
    
    
    // Solve and store in znp1_d, znp1_mp
    
    
    
    
    
    // Exact value of step after Euler prediction
    std::vector<dbl> predicted_d = {dbl(2.19163044775416716032722848776575,-0.203532362690064835307109968099797),
        dbl(1.82976707889226183423896735889061, 2.37621890895572354790101570972854)};
    std::vector<mpfr> predicted_mp = {mpfr("2.19163044775416716032722848776575","-0.203532362690064835307109968099797"),
        mpfr("1.82976707889226183423896735889061", "2.37621890895572354790101570972854")};
    
    // Exact value of step after a two-Newton-step correction
    std::vector<dbl> corrected_d = {dbl(1.23241511572255636681427246524479,0.0106850105579202090198237239792342),
        dbl(0.170574418728320807374606118834491, -0.0023744467906489353377386053287187)};
    std::vector<mpfr> corrected_mp = {mpfr("1.23241511572255636681427246524479","0.0106850105579202090198237239792342"),
        mpfr("0.170574418728320807374606118834491", "-0.0023744467906489353377386053287187")};
    
    
}


BOOST_AUTO_TEST_SUITE_END()

























