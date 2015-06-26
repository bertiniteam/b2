//This file is part of Bertini 2.0.
//
//function_tree_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_tree_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_tree_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  function_tree_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

#include <iostream>

#include <cstdlib>
#include <cmath>
#include <vector>

#include "bertini.hpp"
#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include <eigen3/Eigen/Dense>



using Variable = bertini::Variable;
using Node = bertini::Node;
using Number = bertini::Number;
using SpecialNumber = bertini::SpecialNumber;
using Function = bertini::Function;
using Jacobian = bertini::Jacobian;


int Node::tabcount = 0;


double threshold_clearance_d = 1e-10;

unsigned FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS = 30;
double threshold_clearance_mp = 1e-25;

std::string xstr_real = "3.1";
std::string xstr_imag = "0.3";
std::string ystr_real = "8.8";
std::string ystr_imag = "1.4";
std::string astr_real = "3.4";
std::string astr_imag = "5.6";
std::string bstr_real = "-0.2";
std::string bstr_imag = "-2.1";
std::string zstr_real = "0.8";
std::string zstr_imag = "-1.7";


dbl xnum_dbl(std::stod(xstr_real), std::stod(xstr_imag));
dbl ynum_dbl(std::stod(ystr_real), std::stod(ystr_imag));
dbl anum_dbl(std::stod(astr_real), std::stod(astr_imag));
dbl bnum_dbl(std::stod(bstr_real), std::stod(bstr_imag));
dbl znum_dbl(std::stod(zstr_real), std::stod(zstr_imag));


mpfr xnum_mpfr(xstr_real, xstr_imag);
mpfr ynum_mpfr(ystr_real, ystr_imag);
mpfr anum_mpfr(astr_real, astr_imag);
mpfr bnum_mpfr(bstr_real, bstr_imag);
mpfr znum_mpfr(zstr_real, zstr_imag);

Eigen::Matrix<dbl, 3, 1> var_dbl;
Eigen::Matrix<mpfr, 3, 1> var_mpfr;





BOOST_AUTO_TEST_SUITE(differentiate)

/////////// Basic Operations Alone ///////////////////

BOOST_AUTO_TEST_CASE(just_diff_a_function){
    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    for(auto vv : vars)
    {
        JFunc->EvalJ<dbl>(vv);
        JFunc->EvalJ<mpfr>(vv);
    }
}


BOOST_AUTO_TEST_CASE(diff_3xyz){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = 3*x*y*z;";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {3.0*ynum_dbl*znum_dbl, 3.0*xnum_dbl*znum_dbl, 3.0*ynum_dbl*xnum_dbl};
    std::vector<mpfr> exact_mpfr = {mpfr("3.0")*ynum_mpfr*znum_mpfr,mpfr("3.0")*xnum_mpfr*znum_mpfr,mpfr("3.0")*ynum_mpfr*xnum_mpfr};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_constant){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = 4.5 + i*8.2;";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0, 0.0, 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),mpfr("0.0")};
    
    for(int ii = 0; ii < vars.size(); ++ii)
    {
        BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[ii]).real() - exact_dbl[ii].real() ) < threshold_clearance_d);
        BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[ii]).imag() - exact_dbl[ii].imag()) < threshold_clearance_d);
        BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[ii]).real() - exact_mpfr[ii].real() ) < threshold_clearance_mp);
        BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[ii]).imag() - exact_mpfr[ii].imag() ) < threshold_clearance_mp);
    }
}


BOOST_AUTO_TEST_CASE(diff_sum_xyz_constant){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = x-y+z-4.5+i*7.3;";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {1.0, -1.0, 1.0};
    std::vector<mpfr> exact_mpfr = {mpfr("1.0"),mpfr("-1.0"),mpfr("1.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);

}



BOOST_AUTO_TEST_CASE(diff_x_squared_times_z_cubed){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = (x^2)*(y^3);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {2.0*xnum_dbl*pow(ynum_dbl,3.0), 3.0*pow(ynum_dbl*xnum_dbl,2.0), 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("2.0")*xnum_mpfr*pow(ynum_mpfr,3.0),mpfr("3.0")*pow(ynum_mpfr*xnum_mpfr,2.0),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_x_squared_over_y_cubed){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = (x^2)/(y^3);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {2.0*xnum_dbl/pow(ynum_dbl,3.0), -3.0*pow(xnum_dbl,2.0)/pow(ynum_dbl,4.0), 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("2.0")*xnum_mpfr/pow(ynum_mpfr,3.0),mpfr("-3.0")*pow(xnum_mpfr,2.0)/pow(ynum_mpfr,4.0),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_x_squared_times_lx_plus_numl){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = (x^2)*(x+3);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {3.0*pow(xnum_dbl,2.0) + 6.0*xnum_dbl, 0.0, 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("3.0")*pow(xnum_mpfr,mpfr("2.0")) + mpfr("6.0")*xnum_mpfr,mpfr("0.0"),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}

BOOST_AUTO_TEST_CASE(diff_2y_over_ly_squared_minus_numl){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = y/(y+1);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0, pow(ynum_dbl+1.0,-2.0), 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),pow(ynum_mpfr+mpfr("1.0"),mpfr("-2.0")),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}



BOOST_AUTO_TEST_CASE(diff_sin_x){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = sin(x);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {cos(xnum_dbl), 0.0, 0.0};
    std::vector<mpfr> exact_mpfr = {cos(xnum_mpfr),mpfr("0.0"),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_cos_y){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = cos(y);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0, -1.0*sin(ynum_dbl), 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),-sin(ynum_mpfr),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_tan_z){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = tan(z);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0,0.0, (1.0/cos(znum_dbl))*(1.0/cos(znum_dbl))};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),(mpfr("1.0")/cos(znum_mpfr))*(mpfr("1.0")/cos(znum_mpfr))};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_exp_x){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = exp(x);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {exp(xnum_dbl), 0.0, 0.0};
    std::vector<mpfr> exact_mpfr = {exp(xnum_mpfr),mpfr("0.0"),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_sqrt_y){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = sqrt(y);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0, 0.5/sqrt(ynum_dbl), 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.5")/sqrt(ynum_mpfr),mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}



/////////// Chain Rule ///////////////////
BOOST_AUTO_TEST_CASE(diff_lz_plus_3l_cubed){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = (z+3)^3;";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {0.0, 0.0, 3.0*(pow(znum_dbl+3.0,2.0))};
    std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),mpfr("3.0")*pow(znum_mpfr+mpfr("3.0"),mpfr("2.0"))};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_sin_lx_squared_times_yl){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = sin(x*y);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {cos(xnum_dbl*ynum_dbl)*ynum_dbl, cos(xnum_dbl*ynum_dbl)*xnum_dbl, 0.0};
    std::vector<mpfr> exact_mpfr = {cos(xnum_mpfr*ynum_mpfr)*ynum_mpfr,
        cos(xnum_mpfr*ynum_mpfr)*xnum_mpfr, mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_cos_lx_squaredl){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = cos(x^2);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {-2.0*sin(pow(xnum_dbl,2.0))*xnum_dbl, 0.0, 0.0};
    std::vector<mpfr> exact_mpfr = {mpfr("-2.0")*sin(pow(xnum_mpfr,mpfr("2.0")))*xnum_mpfr,mpfr("0.0"), mpfr("0.0")};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_tan_lx_over_zl){
    using mpfr_float = boost::multiprecision::mpfr_float;
    boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);
    
    std::string str = "function f; variable_group x,y,z; f = tan(x/z);";
    
    bertini::System sys;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    
    var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
    var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
    sys.SetVariables<dbl>(var_dbl);
    sys.SetVariables<mpfr>(var_mpfr);
    
    auto func = sys.function();
    auto vars = sys.variables();
    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
    
    std::vector<dbl> exact_dbl = {1.0/( znum_dbl*pow( cos(xnum_dbl/znum_dbl), 2.0 ) ), 0.0, -xnum_dbl/( pow(znum_dbl, 2.0)*pow( cos(xnum_dbl/znum_dbl), 2.0 ) )};
    std::vector<mpfr> exact_mpfr = {mpfr("1.0")/( znum_mpfr*pow( cos(xnum_mpfr/znum_mpfr), mpfr("2.0") ) ), mpfr("0.0"), -xnum_mpfr/( pow(znum_mpfr, mpfr("2.0"))*pow( cos(xnum_mpfr/znum_mpfr), mpfr("2.0") ) )};
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).real() - exact_dbl[0].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[0]).imag() - exact_dbl[0].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[0]).imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).real() - exact_dbl[1].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[1]).imag() - exact_dbl[1].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[1]).imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);
    
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).real() - exact_dbl[2].real() ) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<dbl>(vars[2]).imag() - exact_dbl[2].imag()) < threshold_clearance_d);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
    BOOST_CHECK(fabs(JFunc->EvalJ<mpfr>(vars[2]).imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}



//BOOST_AUTO_TEST_CASE(diff_sin_lxy_plus_z_squaredl){
//    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
//    
//    bertini::System sys;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    bertini::SystemParser<std::string::const_iterator> S;
//    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
//    auto func = sys.function();
//    auto vars = sys.variables();
//    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
//    for(auto vv : vars)
//    {
//        JFunc->Eval<dbl>(vv);
//        JFunc->Eval<mpfr>(vv);
//    }
//}
//
//
//BOOST_AUTO_TEST_CASE(diff_cos_lxy_plus_z_squaredl){
//    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
//    
//    bertini::System sys;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    bertini::SystemParser<std::string::const_iterator> S;
//    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
//    auto func = sys.function();
//    auto vars = sys.variables();
//    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
//    for(auto vv : vars)
//    {
//        JFunc->Eval<dbl>(vv);
//        JFunc->Eval<mpfr>(vv);
//    }
//}
//
//
//BOOST_AUTO_TEST_CASE(diff_tan_lxy_plus_z_squaredl){
//    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
//    
//    bertini::System sys;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    bertini::SystemParser<std::string::const_iterator> S;
//    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
//    auto func = sys.function();
//    auto vars = sys.variables();
//    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
//    for(auto vv : vars)
//    {
//        JFunc->Eval<dbl>(vv);
//        JFunc->Eval<mpfr>(vv);
//    }
//}
//
//
//BOOST_AUTO_TEST_CASE(diff_exp_lxy_plus_z_squaredl){
//    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
//    
//    bertini::System sys;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    bertini::SystemParser<std::string::const_iterator> S;
//    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
//    auto func = sys.function();
//    auto vars = sys.variables();
//    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
//    for(auto vv : vars)
//    {
//        JFunc->Eval<dbl>(vv);
//        JFunc->Eval<mpfr>(vv);
//    }
//}
//
//
//BOOST_AUTO_TEST_CASE(diff_sqrt_lxy_plus_z_squaredl){
//    std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";
//    
//    bertini::System sys;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    bertini::SystemParser<std::string::const_iterator> S;
//    phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
//    auto func = sys.function();
//    auto vars = sys.variables();
//    auto JFunc = std::make_shared<Jacobian>(func->Differentiate());
//    for(auto vv : vars)
//    {
//        JFunc->Eval<dbl>(vv);
//        JFunc->Eval<mpfr>(vv);
//    }
//}






BOOST_AUTO_TEST_SUITE_END()
