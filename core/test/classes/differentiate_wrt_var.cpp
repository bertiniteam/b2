//This file is part of Bertini 2.
//
//differentiate_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//differentiate_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with differentiate_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Dani Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015, Spring, Summer 2017

#include <iostream>

#include <cstdlib>
#include <cmath>
#include <vector>

#include "bertini2/function_tree.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>




#include "externs.hpp"





BOOST_AUTO_TEST_SUITE(differentiate_wrt_variable)

using bertini::DefaultPrecision;
using dbl = std::complex<double>;
using Variable = bertini::node::Variable;
using Node = bertini::node::Node;
using Function = bertini::node::Function;
using Jacobian = bertini::node::Jacobian;
using bertini::MakeVariable;
using bertini::MakeJacobian;
using dbl = bertini::dbl;
using mpfr = bertini::mpfr;
using mpfr_float = bertini::mpfr_float;
/////////// Basic Operations Alone ///////////////////

BOOST_AUTO_TEST_CASE(just_diff_a_function){
	std::string str = "function f; variable_group x,y,z; f = x*y +y^2 - z*x + 9;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	for(auto vv : vars)
	{
		func->Differentiate(vv);
	}

	std::vector<int> multidegree{1,2,1};
	bool multidegree_ok = multidegree==func->MultiDegree(vars);
	BOOST_CHECK(multidegree_ok);
	
	BOOST_CHECK_EQUAL(func->Degree(vars), 2);
}


BOOST_AUTO_TEST_CASE(diff_3xyz){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = 3*x*y*z;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);


	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();

	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);
	

	std::vector<dbl> exact_dbl = {3.0*ynum_dbl*znum_dbl, 3.0*xnum_dbl*znum_dbl, 3.0*ynum_dbl*xnum_dbl};
	std::vector<mpfr> exact_mpfr = {mpfr("3.0")*ynum_mpfr*znum_mpfr,mpfr("3.0")*xnum_mpfr*znum_mpfr,mpfr("3.0")*ynum_mpfr*xnum_mpfr};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl[0].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl[0].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() / exact_dbl[1].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() / exact_dbl[1].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() / exact_mpfr[1].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() / exact_mpfr[1].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() / exact_dbl[2].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() / exact_dbl[2].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() / exact_mpfr[2].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() / exact_mpfr[2].imag() -1) < threshold_clearance_mp);

	var_mpfr << bertini::complex::rand(),bertini::complex::rand(),bertini::complex::rand();
	sys.SetVariables<mpfr>(var_mpfr);
	exact_mpfr[0] = 3*var_mpfr(1)*var_mpfr(2);
	exact_mpfr[1] = 3*var_mpfr(0)*var_mpfr(2);
	exact_mpfr[2] = 3*var_mpfr(0)*var_mpfr(1);

	dx->Reset();
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	dy->Reset();
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() / exact_mpfr[1].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() / exact_mpfr[1].imag() -1) < threshold_clearance_mp);

	dz->Reset();	
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() / exact_mpfr[2].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() / exact_mpfr[2].imag() -1) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_constant){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = 4.5 + i*8.2;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	std::vector<dbl> exact_dbl = {0.0, 0.0, 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),mpfr("0.0")};

	dx->Reset();
	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real()) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag()) < threshold_clearance_mp);

	dy->Reset();
	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real()) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag()) < threshold_clearance_mp);

	dz->Reset();
	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real()) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag()) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_sum_xyz_constant){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = x-y+z-4.5+i*7.3;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	std::vector<dbl> exact_dbl {1.0, -1.0, 1.0};
	std::vector<mpfr> exact_mpfr {mpfr("1.0"),mpfr("-1.0"),mpfr("1.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);

}



BOOST_AUTO_TEST_CASE(diff_x_squared_times_z_cubed){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = (x^2)*(y^3);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	std::vector<dbl> exact_dbl = {2.0*xnum_dbl*pow(ynum_dbl,3.0), 3.0*pow(ynum_dbl*xnum_dbl,2.0), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("2.0")*xnum_mpfr*pow(ynum_mpfr,3.0),mpfr("3.0")*pow(ynum_mpfr,2)*pow(xnum_mpfr,2.0),mpfr("0.0")};

	

	
	dx->Reset();
	sys.SetVariables<dbl>(var_dbl);
	auto J_val_0_d = dx->Eval<dbl>();

	dy->Reset();
	sys.SetVariables<dbl>(var_dbl);
	auto J_val_1_d = dy->Eval<dbl>();

	dz->Reset();
	sys.SetVariables<dbl>(var_dbl);
	auto J_val_2_d = dz->Eval<dbl>();


	dx->Reset();
	sys.SetVariables<mpfr>(var_mpfr);
	auto J_val_0_mp = dx->Eval<mpfr>();

	dy->Reset();
	sys.SetVariables<mpfr>(var_mpfr);
	auto J_val_1_mp = dy->Eval<mpfr>();

	dz->Reset();
	sys.SetVariables<mpfr>(var_mpfr);
	auto J_val_2_mp = dz->Eval<mpfr>();


	BOOST_CHECK(fabs(J_val_0_d.real() / exact_dbl[0].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_0_d.imag() / exact_dbl[0].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_0_mp.real()/exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(J_val_0_mp.imag()/exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(J_val_1_d.real() / exact_dbl[1].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_1_d.imag() / exact_dbl[1].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_1_mp.real()/exact_mpfr[1].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(J_val_1_mp.imag()/exact_mpfr[1].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(J_val_2_d.real() - exact_dbl[2].real()) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_2_d.imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(J_val_2_mp.real()-exact_mpfr[2].real()) < threshold_clearance_mp);
	BOOST_CHECK(fabs(J_val_2_mp.imag()-exact_mpfr[2].imag()) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_x_squared_over_y_cubed){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = (x^2)/(y^3);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);
	



	std::vector<dbl> exact_dbl = {2.0*xnum_dbl/pow(ynum_dbl,3.0), -3.0*pow(xnum_dbl,2.0)/pow(ynum_dbl,4.0), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("2.0")*xnum_mpfr/pow(ynum_mpfr,3.0),mpfr("-3.0")*pow(xnum_mpfr,2.0)/pow(ynum_mpfr,4.0),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_x_squared_times_lx_plus_numl){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = (x^2)*(x+3);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);





	std::vector<dbl> exact_dbl = {3.0*pow(xnum_dbl,2.0) + 6.0*xnum_dbl, 0.0, 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("3.0")*pow(xnum_mpfr,mpfr("2.0")) + mpfr("6.0")*xnum_mpfr,mpfr("0.0"),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl[0].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl[0].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}

BOOST_AUTO_TEST_CASE(diff_2y_over_ly_squared_minus_numl){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = y/(y+1);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	

	std::vector<dbl> exact_dbl = {0.0, pow(ynum_dbl+1.0,-2.0), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),pow(ynum_mpfr+mpfr("1.0"),mpfr("-2.0")),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}



BOOST_AUTO_TEST_CASE(diff_sin_x){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = sin(x);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);



	
	std::vector<dbl> exact_dbl = {cos(xnum_dbl), 0.0, 0.0};
	std::vector<mpfr> exact_mpfr = {cos(xnum_mpfr),mpfr("0.0"),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_cos_y){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = cos(y);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	

	std::vector<dbl> exact_dbl = {0.0, -1.0*sin(ynum_dbl), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),-sin(ynum_mpfr),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_tan_z){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = tan(z);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);



	std::vector<dbl> exact_dbl = {0.0,0.0, (1.0/cos(znum_dbl))*(1.0/cos(znum_dbl))};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),(mpfr("1.0")/cos(znum_mpfr))*(mpfr("1.0")/cos(znum_mpfr))};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_exp_x){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = exp(x);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	std::vector<dbl> exact_dbl = {exp(xnum_dbl), 0.0, 0.0};
	std::vector<mpfr> exact_mpfr = {exp(xnum_mpfr),mpfr("0.0"),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_log_x){
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = log(x^2+y);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);
	sys.Differentiate();

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);


	std::vector<dbl> exact_dbl = {2.0*xnum_dbl/(xnum_dbl*xnum_dbl+ynum_dbl), 1.0/(xnum_dbl*xnum_dbl+ynum_dbl), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr(2.0)*xnum_mpfr/(xnum_mpfr*xnum_mpfr+ynum_mpfr), mpfr(1.0)/(xnum_mpfr*xnum_mpfr+ynum_mpfr),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_sqrt_y){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = sqrt(y);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	

	std::vector<dbl> exact_dbl = {0.0, 0.5/sqrt(ynum_dbl), 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.5")/sqrt(ynum_mpfr),mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}



/////////// Chain Rule ///////////////////
BOOST_AUTO_TEST_CASE(diff_lz_plus_3l_cubed){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = (z+3)^3;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();

	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	std::vector<dbl> exact_dbl = {0.0, 0.0, 3.0*(pow(znum_dbl+3.0,2.0))};
	std::vector<mpfr> exact_mpfr = {mpfr("0.0"),mpfr("0.0"),mpfr("3.0")*pow(znum_mpfr+mpfr("3.0"),mpfr("2.0"))};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}




BOOST_AUTO_TEST_CASE(diff_x_squared_plus_y_squared_plus_z_squared){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = x^2+y^2+z^2;";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	


	std::vector<dbl> exact_dbl = {2.0*xnum_dbl, 2.0*ynum_dbl, 2.0*znum_dbl};
	std::vector<mpfr> exact_mpfr = {mpfr("2.0")*xnum_mpfr, mpfr("2.0")*ynum_mpfr, mpfr("2.0")*znum_mpfr};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() / exact_dbl[1].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() / exact_dbl[1].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}





BOOST_AUTO_TEST_CASE(diff_sin_lx_squared_times_yl){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = sin(x*y);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	std::vector<dbl> exact_dbl = {cos(xnum_dbl*ynum_dbl)*ynum_dbl, cos(xnum_dbl*ynum_dbl)*xnum_dbl, 0.0};
	std::vector<mpfr> exact_mpfr = {cos(xnum_mpfr*ynum_mpfr)*ynum_mpfr,
		cos(xnum_mpfr*ynum_mpfr)*xnum_mpfr, mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl[0].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl[0].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() / exact_dbl[1].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() / exact_dbl[1].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() / exact_mpfr[1].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() / exact_mpfr[1].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_cos_lx_squaredl){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = cos(x^2);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	std::vector<dbl> exact_dbl = {-2.0*sin(pow(xnum_dbl,2.0))*xnum_dbl, 0.0, 0.0};
	std::vector<mpfr> exact_mpfr = {mpfr("-2.0")*sin(pow(xnum_mpfr,mpfr("2.0")))*xnum_mpfr,mpfr("0.0"), mpfr("0.0")};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl[0].real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl[0].imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr[0].real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr[0].imag() -1) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(diff_tan_lx_over_zl){
    bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x,y,z; f = tan(x/z);";

	bertini::System sys;
	bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Eigen::Matrix<dbl, 3, 1> var_dbl;
	Eigen::Matrix<mpfr, 3, 1> var_mpfr;

	var_dbl << xnum_dbl, ynum_dbl, znum_dbl;
	var_mpfr << xnum_mpfr, ynum_mpfr, znum_mpfr;
	sys.SetVariables<dbl>(var_dbl);
	sys.SetVariables<mpfr>(var_mpfr);

	auto func = sys.Function(0);
	auto vars = sys.Variables();
	
	auto dx = func->Differentiate(vars[0]);
	auto dy = func->Differentiate(vars[1]);
	auto dz = func->Differentiate(vars[2]);

	std::vector<dbl> exact_dbl = {1.0/( znum_dbl*pow( cos(xnum_dbl/znum_dbl), 2.0 ) ), 0.0, -xnum_dbl/( pow(znum_dbl, 2.0)*pow( cos(xnum_dbl/znum_dbl), 2.0 ) )};
	std::vector<mpfr> exact_mpfr = {mpfr("1.0")/( znum_mpfr*pow( cos(xnum_mpfr/znum_mpfr), mpfr("2.0") ) ), mpfr("0.0"), -xnum_mpfr/( pow(znum_mpfr, mpfr("2.0"))*pow( cos(xnum_mpfr/znum_mpfr), mpfr("2.0") ) )};

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() - exact_dbl[0].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() - exact_dbl[0].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() - exact_mpfr[0].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() - exact_mpfr[0].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dy->Eval<dbl>().real() - exact_dbl[1].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<dbl>().imag() - exact_dbl[1].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().real() - exact_mpfr[1].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dy->Eval<mpfr>().imag() - exact_mpfr[1].imag() ) < threshold_clearance_mp);

	BOOST_CHECK(fabs(dz->Eval<dbl>().real() - exact_dbl[2].real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<dbl>().imag() - exact_dbl[2].imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().real() - exact_mpfr[2].real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dz->Eval<mpfr>().imag() - exact_mpfr[2].imag() ) < threshold_clearance_mp);
}





BOOST_AUTO_TEST_CASE(arcsine_differentiate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = asin(pow(x,2)+1);
	auto dx = N->Differentiate(x);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	//(2*x)/(1 - (x^2 + 1)^2)^(1/2)
	dbl exact_dbl = 2.0*xnum_dbl / pow(1.0 - pow((xnum_dbl*xnum_dbl + 1.0),2),0.5);
	mpfr exact_mpfr = bertini::complex(2.0)*xnum_mpfr / pow(bertini::complex(1.0) - pow(xnum_mpfr*xnum_mpfr + bertini::complex(1.0),2),mpfr(0.5));

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(arccosine_differentiate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = acos(pow(x,2)+1);
	auto dx = N->Differentiate(x);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	dbl exact_dbl = -2.0*xnum_dbl / pow(1.0 - pow((xnum_dbl*xnum_dbl + 1.0),2),0.5);
	mpfr exact_mpfr = -bertini::complex(2.0)*xnum_mpfr / pow(bertini::complex(1.0) - pow(xnum_mpfr*xnum_mpfr + bertini::complex(1.0),2),mpfr(0.5));

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
}

BOOST_AUTO_TEST_CASE(arctangent_differentiate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = atan(pow(x,2)+1);
	auto dx = N->Differentiate(x);


	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	//(2*x)/((x^2 + 1)^2 + 1)
	dbl exact_dbl = 2.0*xnum_dbl / ( pow((xnum_dbl*xnum_dbl + 1.0),2) + 1.0);
	mpfr exact_mpfr = bertini::complex(2.0)*xnum_mpfr / ( pow((xnum_mpfr*xnum_mpfr + bertini::complex(1.0)),2) + bertini::complex(1.0));

	BOOST_CHECK(fabs(dx->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(dx->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(integer_power)
{
	

	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> t = MakeVariable("t");

	auto f = pow(x - 1,2)*(1-t) + (pow(x,2) + 1)*t;
	auto dx = f->Differentiate(x);


	x->set_current_value(mpfr("-0.844487","-0.535576"));
	t->set_current_value(mpfr("0.779871","0.712645"));

	auto J = dx->Eval<mpfr>();
	BOOST_CHECK(abs(J - mpfr("-2.129232","0.354138")) < threshold_clearance_mp);


	x->set_current_value(mpfr("0.900000000000000","0.435889894354067355223698198386"));
    t->set_current_value(mpfr("0.1"));
    dx->Reset();
	J = dx->Eval<mpfr>();

	BOOST_CHECK(abs( real(J) ) < threshold_clearance_mp);
	BOOST_CHECK(abs( imag(J) - mpfr_float("0.871779788708134710447396396772")) < threshold_clearance_mp);
}



BOOST_AUTO_TEST_SUITE_END()
