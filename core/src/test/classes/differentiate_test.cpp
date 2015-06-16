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

#include "bertini.hpp"
#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>



using Variable = bertini::Variable;
using Node = bertini::Node;
using Number = bertini::Number;
using SpecialNumber = bertini::SpecialNumber;
using Function = bertini::Function;


int Node::tabcount = 0;


double threshold_clearance_d = 1e-15;

unsigned FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS = 30;
double threshold_clearance_mp = 1e-27;

std::string xstr_real = "3.1";
std::string xstr_imag = "4.1";
std::string ystr_real = "8.8";
std::string ystr_imag = "9.9";
std::string astr_real = "3.4";
std::string astr_imag = "5.6";
std::string bstr_real = "-0.2";
std::string bstr_imag = "-2.1";
std::string pstr_real = "0.8";
std::string pstr_imag = "-1.7";


dbl xnum_dbl(std::stod(xstr_real), std::stod(xstr_imag));
dbl ynum_dbl(std::stod(ystr_real), std::stod(ystr_imag));
dbl anum_dbl(std::stod(astr_real), std::stod(astr_imag));
dbl bnum_dbl(std::stod(bstr_real), std::stod(bstr_imag));
dbl pnum_dbl(std::stod(pstr_real), std::stod(pstr_imag));

mpfr xnum_mpfr(xstr_real, xstr_imag);
mpfr ynum_mpfr(ystr_real, ystr_imag);
mpfr anum_mpfr(astr_real, astr_imag);
mpfr bnum_mpfr(bstr_real, bstr_imag);
mpfr pnum_mpfr(pstr_real, pstr_imag);





BOOST_AUTO_TEST_SUITE(differentiate)

/////////// Basic Operations Alone ///////////////////

BOOST_AUTO_TEST_CASE(manual_construction_num_squared){
    std::string str = "function f; variable_group x1, x2; f = -x1;";
    
    bertini::System sys;
    std::shared_ptr<Function> func, diffFunc;
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bertini::SystemParser<std::string::const_iterator> S;
    bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
    func = sys.function();
    diffFunc = std::make_shared<Function>(func->Differentiate());
    diffFunc->print(std::cout);
}







BOOST_AUTO_TEST_SUITE_END()
