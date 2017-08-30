//This file is part of Bertini 2.
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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

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

#include "bertini2/function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include "externs.hpp"



BOOST_AUTO_TEST_SUITE(function_tree_class)


using mpq_rational = bertini::mpq_rational;

using Variable = bertini::node::Variable;
using Node = bertini::node::Node;
using Float = bertini::node::Float;

using dbl = bertini::dbl;
using mpfr = bertini::mpfr;

using bertini::MakeVariable;
using bertini::MakeFloat;

using namespace bertini;

/////////// Basic Operations Alone ///////////////////

BOOST_AUTO_TEST_CASE(manual_construction_num_squared){
	using mpfr_float = bertini::mpfr_float;
	
	
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = anum_dbl*anum_dbl;
	mpfr exact_mpfr = anum_mpfr*anum_mpfr;
	
	std::shared_ptr<Node> N = a;
	
	N *= N;
	BOOST_CHECK_EQUAL(N->Degree(),0);
	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_x_squared){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	dbl exact_dbl = xnum_dbl*xnum_dbl;
	mpfr exact_mpfr = xnum_mpfr*xnum_mpfr;
	
	std::shared_ptr<Node> N = x;
	
	N *= N;
	
	x->set_current_value<dbl>(std::complex<double>(xnum_dbl));
	x->set_current_value<mpfr>(bertini::complex(xnum_mpfr));
	
	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),2);
	

}



BOOST_AUTO_TEST_CASE(default_constructed_variable_is_not_nan){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	using bertini::isnan;
	BOOST_CHECK(!isnan(x->Eval<dbl>()));
	BOOST_CHECK(!isnan(x->Eval<mpfr>()));
}

BOOST_AUTO_TEST_CASE(default_constructed_variable_is_not_zero){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	BOOST_CHECK(x->Eval<dbl>()!=dbl(0,0));
	BOOST_CHECK(x->Eval<mpfr>()!=mpfr(0,0));
}


BOOST_AUTO_TEST_CASE(default_constructed_variable_is_not_one){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	BOOST_CHECK(x->Eval<dbl>()!=dbl(1,0));
	BOOST_CHECK(x->Eval<mpfr>()!=mpfr(1,0));
}

BOOST_AUTO_TEST_CASE(self_multiplication){

	auto rat_coeff = [](){return bertini::MakeRational(bertini::node::Rational::Rand());};

	std::shared_ptr<Node> v = rat_coeff();
	std::shared_ptr<Node> N = v*v;
}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_x){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	dbl exact_dbl = sqrt(xnum_dbl);
	mpfr exact_mpfr = sqrt(xnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x, mpq_rational(1,2));
	
	x->set_current_value<dbl>(std::complex<double>(xnum_dbl));
	x->set_current_value<mpfr>(bertini::complex(xnum_mpfr));
	
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_plus_y_plus_number){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+ynum_dbl+anum_dbl;
	mpfr exact_mpfr = xnum_mpfr+ynum_mpfr+anum_mpfr;
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = x+y+a;
	
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = a+x+y;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = y+a+x;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	
	N = y+x+a;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	
	N = x+a+y;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = a+y+x;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_x_minus_y_minus_number){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl-ynum_dbl-anum_dbl;
	mpfr exact_mpfr = xnum_mpfr-ynum_mpfr-anum_mpfr;
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = x-y-a;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	N = x-a-y;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);


	N = -a-y+x;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);


	N = -a+x-y;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	N = -y-a+x;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);


	N = -y+x-a;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);


	
}


BOOST_AUTO_TEST_CASE(manual_construction_x_times_y_times_number){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl*ynum_dbl*anum_dbl;
	mpfr exact_mpfr = xnum_mpfr*ynum_mpfr*anum_mpfr;
	
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = x*y*a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);
	
	BOOST_CHECK(N->IsPolynomial());
	
	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	N = a*x*y;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	N = y*a*x;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	
	N = y*x*a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	
	N = x*a*y;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
	N = a*y*x;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_x_divide_y){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl/ynum_dbl;
	mpfr exact_mpfr = xnum_mpfr/ynum_mpfr;
	
	std::shared_ptr<Node> N = x/y;
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	exact_dbl = ynum_dbl/xnum_dbl;
	exact_mpfr = ynum_mpfr/xnum_mpfr;
	
	N = y/x;

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_negate_x){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = -xnum_dbl;
	mpfr exact_mpfr = -xnum_mpfr;
	
	std::shared_ptr<Node> N = -x;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}







/////////// Basic Operations Combined ///////////////////


BOOST_AUTO_TEST_CASE(manual_construction_lx_plus_y_plus_num1l_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl+ynum_dbl+anum_dbl,pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr+ynum_mpfr+anum_mpfr,pnum_mpfr);
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = pow(x+y+a,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);

}


BOOST_AUTO_TEST_CASE(manual_construction_lx_minus_y_minus_num1l_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl-ynum_dbl-anum_dbl,pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr-ynum_mpfr-anum_mpfr,pnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x-y-a,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_times_y_times_num1l_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl*ynum_dbl*anum_dbl,pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr*ynum_mpfr*anum_mpfr,pnum_mpfr);
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = pow(x*y*a,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag()) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = pow(x,p)*pow(y,p)*pow(a,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	exact_dbl = pow(xnum_dbl,pnum_dbl)*pow(ynum_dbl,pnum_dbl)*pow(anum_dbl,pnum_dbl);
	exact_mpfr = pow(xnum_mpfr,pnum_mpfr)*pow(ynum_mpfr,pnum_mpfr)*pow(anum_mpfr,pnum_mpfr);

	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_over_yl_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl/ynum_dbl,pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr/ynum_mpfr,pnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x/y,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = pow(x,p)/pow(y,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_lnegative_xl_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(-xnum_dbl,pnum_dbl);
	mpfr exact_mpfr = pow(-xnum_mpfr,pnum_mpfr);
	
	std::shared_ptr<Node> N = -x;
	N = pow(N,p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_negate_x_plus_y_plus_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl+ynum_dbl+anum_dbl);
	mpfr exact_mpfr = -(xnum_mpfr+ynum_mpfr+anum_mpfr);
	
	std::shared_ptr<Node> N = -(x+y+a);
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = -x-y-a;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_negate_x_minus_y_minus_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl-ynum_dbl-anum_dbl);
	mpfr exact_mpfr = -(xnum_mpfr-ynum_mpfr-anum_mpfr);
	
	std::shared_ptr<Node> N = -(x-y-a);
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = -x+y+a;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}



BOOST_AUTO_TEST_CASE(manual_construction_negate_x_times_y_times_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl*ynum_dbl*anum_dbl);
	mpfr exact_mpfr = -(xnum_mpfr*ynum_mpfr*anum_mpfr);
	
	std::shared_ptr<Node> N = -(x*y*a);
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = (-x)*y*a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = x*(-y)*a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = x*y*(-a);
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}



BOOST_AUTO_TEST_CASE(manual_construction_negate_x_over_y){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl/ynum_dbl);
	mpfr exact_mpfr = -(xnum_mpfr/ynum_mpfr);
	
	std::shared_ptr<Node> N = -(x/y);
	BOOST_CHECK_EQUAL(N->Degree(),-1); 
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	
	N = -x/y;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);


	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	
	N = x/(-y);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	;
	
}


BOOST_AUTO_TEST_CASE(manual_construction_negate_x_pow_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = -pow(xnum_dbl,pnum_dbl);
	mpfr exact_mpfr = -pow(xnum_mpfr,pnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x,p);
	N = -N;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	
	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}







/////////// Order of Operations ///////////////////

BOOST_AUTO_TEST_CASE(manual_construction_x_times_y_over_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl*ynum_dbl/anum_dbl;
	mpfr exact_mpfr = xnum_mpfr*ynum_mpfr/anum_mpfr;
	
	std::shared_ptr<Node> N = x*y/a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	


	N = (x*y)/a;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	N = x/a*y;
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK_SMALL(N->Eval<dbl>().real()/exact_dbl.real()-1.0,  threshold_clearance_d);
	BOOST_CHECK_SMALL(N->Eval<dbl>().imag()/exact_dbl.imag()-1.0,  threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_plus_num1l_times_ly_plus_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = (xnum_dbl+anum_dbl)*(ynum_dbl+bnum_dbl);
	mpfr exact_mpfr = (xnum_mpfr+anum_mpfr)*(ynum_mpfr+bnum_mpfr);
	
	std::shared_ptr<Node> N = (x+a)*(y+b);
	BOOST_CHECK_EQUAL(N->Degree(),2);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(manual_construction_x_plus_num1_times_y_plus_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+anum_dbl*ynum_dbl+bnum_dbl;
	mpfr exact_mpfr = xnum_mpfr+anum_mpfr*ynum_mpfr+bnum_mpfr;
	
	std::shared_ptr<Node> N = x+a*y+b;
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_plus_num1l_over_ly_plus_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = (xnum_dbl+anum_dbl)/(ynum_dbl+bnum_dbl);
	mpfr exact_mpfr = (xnum_mpfr+anum_mpfr)/(ynum_mpfr+bnum_mpfr);
	
	std::shared_ptr<Node> N = (x+a)/(y+b);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_plus_num1_over_y_plus_num2){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Variable> y = MakeVariable("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+anum_dbl/ynum_dbl+bnum_dbl;
	mpfr exact_mpfr = xnum_mpfr+anum_mpfr/ynum_mpfr+bnum_mpfr;
	
	std::shared_ptr<Node> N = x+a/y+b;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);
	BOOST_CHECK_EQUAL(N->Degree(y),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_pow_num2l_plus_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)+anum_dbl;
	mpfr exact_mpfr = pow(xnum_mpfr,pnum_mpfr)+anum_mpfr;
	
	std::shared_ptr<Node> N = pow(x,p)+a;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_plus_lnum1_pow_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(anum_dbl,pnum_dbl)+xnum_dbl;
	mpfr exact_mpfr = pow(anum_mpfr,pnum_mpfr)+xnum_mpfr;
	
	std::shared_ptr<Node> N = x+pow(a,p);
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_times_lnum1_pow_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(anum_dbl,pnum_dbl)*xnum_dbl;
	mpfr exact_mpfr = pow(anum_mpfr,pnum_mpfr)*xnum_mpfr;
	
	std::shared_ptr<Node> N = x*pow(a,p);
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_pow_num2l_times_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)*anum_dbl;
	mpfr exact_mpfr = pow(xnum_mpfr,pnum_mpfr)*anum_mpfr;
	
	std::shared_ptr<Node> N = pow(x,p)*a;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_lx_pow_num2l_over_num1){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)/anum_dbl;
	mpfr exact_mpfr = pow(xnum_mpfr,pnum_mpfr)/anum_mpfr;
	
	std::shared_ptr<Node> N = pow(x,p)/a;
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(! N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_pow_lsqrt_xl_num)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	auto exact_dbl = pow(sqrt(xnum_dbl),anum_dbl);
	auto exact_mpfr = pow(sqrt(xnum_mpfr),anum_mpfr);

	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(xnum_mpfr);

	std::shared_ptr<Node> N = pow(sqrt(x),a);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	N = pow(x,mpq_rational(1,2)); // N = sqrt(x)
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	N = pow(N,a); // N = sqrt(x)^a
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK_SMALL(N->Eval<dbl>().real() / exact_dbl.real() -1, threshold_clearance_d);
	BOOST_CHECK_SMALL(N->Eval<dbl>().imag() - exact_dbl.imag(), 10*threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

}


BOOST_AUTO_TEST_CASE(manual_construction_x_over_lnum1_pow_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = xnum_dbl/pow(anum_dbl,pnum_dbl);
	mpfr exact_mpfr = xnum_mpfr/pow(anum_mpfr,pnum_mpfr);
	
	std::shared_ptr<Node> N = x/pow(a,p);
	BOOST_CHECK_EQUAL(N->Degree(),1);
	BOOST_CHECK_EQUAL(N->Degree(x),1);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_pow_lnum1_plus_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,anum_dbl+pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr,anum_mpfr+pnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x,a+p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_pow_lnum1_times_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,anum_dbl*pnum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr,anum_mpfr*pnum_mpfr);
	
	std::shared_ptr<Node> N = pow(x,a*p);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_pow_lnum1_over_num2l){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl/anum_dbl);
	mpfr exact_mpfr = pow(xnum_mpfr,pnum_mpfr/anum_mpfr);
	
	std::shared_ptr<Node> N = pow(x,p/a);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);


	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}









///////////// Special Functions ///////////////////
BOOST_AUTO_TEST_CASE(manual_construction_sin_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = sin(anum_dbl);
	mpfr exact_mpfr = sin(anum_mpfr);
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = sin(a);
	
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);

}


BOOST_AUTO_TEST_CASE(manual_construction_cos_num){
	using mpfr_float = bertini::mpfr_float;
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = cos(anum_dbl);
	mpfr exact_mpfr = cos(anum_mpfr);
	
	std::shared_ptr<Node> N = cos(a);
	
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_tan_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = tan(anum_dbl);
	mpfr exact_mpfr = tan(anum_mpfr);
	
	std::shared_ptr<Node> N = tan(a);
	
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_exp_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = exp(anum_dbl);
	mpfr exact_mpfr = exp(anum_mpfr);
	
	dbl temp_d;
	mpfr temp_mp;
	
	std::shared_ptr<Node> N = exp(a);
	
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N->EvalInPlace<dbl>(temp_d);
	N->EvalInPlace<mpfr>(temp_mp);
	BOOST_CHECK(fabs(temp_d.real() / exact_dbl.real() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_d.imag() / exact_dbl.imag() -1) < threshold_clearance_d);
	BOOST_CHECK(fabs(temp_mp.real() / exact_mpfr.real() -1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(temp_mp.imag() / exact_mpfr.imag() -1) < threshold_clearance_mp);

}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = sqrt(anum_dbl);
	mpfr exact_mpfr = sqrt(anum_mpfr);
	
	std::shared_ptr<Node> N = sqrt(a);
	
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_sin_of_lx_plus_numl){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = sin(xnum_dbl+anum_dbl);
	mpfr exact_mpfr = sin(xnum_mpfr+anum_mpfr);
	
	std::shared_ptr<Node> N = sin(x+a);
	
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_cos_of_lx_times_numl){
	using mpfr_float = bertini::mpfr_float;
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = cos(xnum_dbl*anum_dbl);
	mpfr exact_mpfr = cos(xnum_mpfr*anum_mpfr);
	
	std::shared_ptr<Node> N = cos(x*a);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_tan_of_lx_over_numl){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = tan(xnum_dbl/anum_dbl);
	mpfr exact_mpfr = tan(xnum_mpfr/anum_mpfr);
	
	std::shared_ptr<Node> N = tan(x/a);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());
	
	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_exp_of_negative_num){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = exp(-anum_dbl);
	mpfr exact_mpfr = exp(-anum_mpfr);
	
	std::shared_ptr<Node> N = exp(-a);
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());
	
	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_of_lx_pow_numl){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = sqrt(pow(xnum_dbl,anum_dbl));
	mpfr exact_mpfr = sqrt(pow(xnum_mpfr,anum_mpfr));
	
	std::shared_ptr<Node> N = sqrt(pow(x,a));
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(!N->IsPolynomial());
	
	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
	
	N = pow(x,a);
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);
	N = pow(N,mpq_rational(1,2));
	BOOST_CHECK_EQUAL(N->Degree(),-1);
	BOOST_CHECK_EQUAL(N->Degree(x),-1);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);
}






BOOST_AUTO_TEST_CASE(arcsine_evaluate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = asin(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	dbl exact_dbl = asin(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr = asin(pow(xnum_mpfr,2)+bertini::complex(1.0));

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(arccosine_evaluate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = acos(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	dbl exact_dbl = acos(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr = acos(pow(xnum_mpfr,2)+bertini::complex(1.0));

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(arctangent_evaluate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = atan(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	dbl exact_dbl = atan(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr = atan(pow(xnum_mpfr,2)+bertini::complex(1.0));

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(log_evaluate)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	auto N = log(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::complex(xstr_real,xstr_imag));

	dbl exact_dbl = log(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr = log(pow(xnum_mpfr,2)+bertini::complex(1.0));

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}


///////////// Special Numbers ///////////////////
BOOST_AUTO_TEST_CASE(manual_construction_pi){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	using std::atan;
	dbl exact_dbl(4*atan(1.0),0);
	mpfr exact_mpfr = mpfr_float("4.0")*atan(mpfr_float("1.0"));
	
	auto N = bertini::node::Pi();
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK_EQUAL(N->Eval<dbl>().imag(),0.0);

	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK_EQUAL(N->Eval<mpfr>().imag(),0.0);

	BOOST_CHECK(N->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(manual_construction_e){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	dbl exact_dbl(exp(1.0),0);
	mpfr exact_mpfr = exp(mpfr_float("1.0"));
	
	auto N = bertini::node::E();
	BOOST_CHECK_EQUAL(N->Degree(),0);
	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(N->Eval<dbl>().imag() == 0.0);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(N->Eval<mpfr>().imag() - exact_mpfr.imag() == 0.0);
}


BOOST_AUTO_TEST_CASE(manual_construction_i){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	dbl exact_dbl(0.0,1.0);
	mpfr exact_mpfr = mpfr("0.0","1.0");
	
	auto N = bertini::node::I();
	BOOST_CHECK_EQUAL(N->Degree(),0);

	BOOST_CHECK(N->IsPolynomial());

	BOOST_CHECK(N->Eval<dbl>().real() == 0.0);
	BOOST_CHECK(N->Eval<dbl>().imag() == 1.0);
	BOOST_CHECK(N->Eval<mpfr>().real() == 0.0);
	BOOST_CHECK(N->Eval<mpfr>().imag() == 1.0);
}





BOOST_AUTO_TEST_CASE(function_tree_combine_product_of_two_integer_powers)
{
	std::shared_ptr<Variable> x = MakeVariable("x");
	std::shared_ptr<Node> N,M,P;

	N = pow(x,5);
	M = pow(x,2);

	P = N*M;
	BOOST_CHECK_EQUAL(P->Degree(), 7);

	BOOST_CHECK(std::dynamic_pointer_cast<bertini::node::IntegerPowerOperator>(P));

	P = N/M;
	BOOST_CHECK(!std::dynamic_pointer_cast<bertini::node::IntegerPowerOperator>(P));
	BOOST_CHECK_EQUAL(P->Degree(), -1);

	N = pow(x,3);
	M = pow(x,-1);

	P = N/M;
	BOOST_CHECK(!std::dynamic_pointer_cast<bertini::node::IntegerPowerOperator>(P));
	BOOST_CHECK_EQUAL(P->Degree(), -1);
}

BOOST_AUTO_TEST_CASE(long_arithmetic_chain)
{

	std::vector<std::shared_ptr<Node>> polytypes(9);

	auto rat_coeff = [](){return bertini::MakeRational(bertini::node::Rational::Rand());};
	for (unsigned int ii=0; ii<9; ++ii)
	    polytypes[ii] = rat_coeff();

	std::shared_ptr<Node> v = rat_coeff();
	std::shared_ptr<Node> m12_2 = rat_coeff();
	std::shared_ptr<Node> cosb = rat_coeff();
	std::shared_ptr<Node> sinb = rat_coeff();
	std::shared_ptr<Node> sin2b = rat_coeff();

	std::shared_ptr<Node> N = 2*m12_2/sin2b - v*v * polytypes[9-1] / 16;

}

BOOST_AUTO_TEST_SUITE_END()






