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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.

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
using mpfr = bertini::mpfr_complex;

using bertini::Variable::Make;
using bertini::MakeFloat;

template<typename NumType> using Mat = bertini::Mat<NumType>;

using namespace bertini;

/////////// Basic Operations Alone ///////////////////

BOOST_AUTO_TEST_CASE(manual_construction_num_squared){
	using mpfr_float = bertini::mpfr_float;
	
	
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	dbl exact_dbl = anum_dbl*anum_dbl;
	mpfr exact_mpfr{anum_mpfr*anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	dbl exact_dbl = xnum_dbl*xnum_dbl;
	mpfr exact_mpfr{xnum_mpfr*xnum_mpfr};
	
	std::shared_ptr<Node> N = x;
	
	N *= N;
	
	x->set_current_value<dbl>(std::complex<double>(xnum_dbl));
	x->set_current_value<mpfr>(bertini::mpfr_complex(xnum_mpfr));
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	using std::isnan;
	using bertini::isnan;
	BOOST_CHECK(!isnan(x->Eval<dbl>()));
	BOOST_CHECK(!isnan(x->Eval<mpfr>()));
}

BOOST_AUTO_TEST_CASE(default_constructed_variable_is_not_zero){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	BOOST_CHECK(x->Eval<dbl>()!=dbl(0,0));
	BOOST_CHECK(x->Eval<mpfr>()!=mpfr(0,0));
}


BOOST_AUTO_TEST_CASE(default_constructed_variable_is_not_one){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	BOOST_CHECK(x->Eval<dbl>()!=dbl(1,0));
	BOOST_CHECK(x->Eval<mpfr>()!=mpfr(1,0));
}

BOOST_AUTO_TEST_CASE(self_multiplication){

	auto rat_coeff = [](){return bertini::Rational::Make(bertini::node::Rational::Rand());};

	std::shared_ptr<Node> v = rat_coeff();
	std::shared_ptr<Node> N = v*v;
}


BOOST_AUTO_TEST_CASE(rational_node_eval_sane_precision_random_rat){

	DefaultPrecision(16);

	std::shared_ptr<Node> frac = bertini::Rational::Make(bertini::node::Rational::Rand());
	mpfr_complex result = frac->Eval<mpfr_complex>();

	BOOST_CHECK_EQUAL(Precision(result),16);
}

BOOST_AUTO_TEST_CASE(rational_node_eval_sane_precision_one_half){

	DefaultPrecision(16);

	std::shared_ptr<Node> frac = bertini::Rational::Make(mpq_rational(1,2));
	mpfr_complex result = frac->Eval<mpfr_complex>();

	BOOST_CHECK_EQUAL(Precision(result),16);
}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_x){
	using mpfr_float = bertini::mpfr_float;
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	dbl exact_dbl = sqrt(xnum_dbl);
	mpfr exact_mpfr{sqrt(xnum_mpfr)};
	
	std::shared_ptr<Node> N = pow(x, mpq_rational(1,2));
	
	x->set_current_value<dbl>(std::complex<double>(xnum_dbl));
	x->set_current_value<mpfr>(bertini::mpfr_complex(xnum_mpfr));
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+ynum_dbl+anum_dbl;
	mpfr exact_mpfr{xnum_mpfr+ynum_mpfr+anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl-ynum_dbl-anum_dbl;
	mpfr exact_mpfr{xnum_mpfr-ynum_mpfr-anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl*ynum_dbl*anum_dbl;
	mpfr exact_mpfr{xnum_mpfr*ynum_mpfr*anum_mpfr};
	
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl/ynum_dbl;
	mpfr exact_mpfr{xnum_mpfr/ynum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = -xnum_dbl;
	mpfr exact_mpfr{-xnum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl+ynum_dbl+anum_dbl,pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr+ynum_mpfr+anum_mpfr,pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl-ynum_dbl-anum_dbl,pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr-ynum_mpfr-anum_mpfr,pnum_mpfr)};
	
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
	
	
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl*ynum_dbl*anum_dbl,pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr*ynum_mpfr*anum_mpfr,pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = pow(xnum_dbl/ynum_dbl,pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr/ynum_mpfr,pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(-xnum_dbl,pnum_dbl);
	mpfr exact_mpfr{pow(-xnum_mpfr,pnum_mpfr)};
	
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
	
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl+ynum_dbl+anum_dbl);
	mpfr exact_mpfr{-(xnum_mpfr+ynum_mpfr+anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl-ynum_dbl-anum_dbl);
	mpfr exact_mpfr{-(xnum_mpfr-ynum_mpfr-anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl*ynum_dbl*anum_dbl);
	mpfr exact_mpfr{-(xnum_mpfr*ynum_mpfr*anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = -(xnum_dbl/ynum_dbl);
	mpfr exact_mpfr{-(xnum_mpfr/ynum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = -pow(xnum_dbl,pnum_dbl);
	mpfr exact_mpfr{-pow(xnum_mpfr,pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl*ynum_dbl/anum_dbl;
	mpfr exact_mpfr{xnum_mpfr*ynum_mpfr/anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = (xnum_dbl+anum_dbl)*(ynum_dbl+bnum_dbl);
	mpfr exact_mpfr{(xnum_mpfr+anum_mpfr)*(ynum_mpfr+bnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+anum_dbl*ynum_dbl+bnum_dbl;
	mpfr exact_mpfr{xnum_mpfr+anum_mpfr*ynum_mpfr+bnum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = (xnum_dbl+anum_dbl)/(ynum_dbl+bnum_dbl);
	mpfr exact_mpfr{(xnum_mpfr+anum_mpfr)/(ynum_mpfr+bnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Variable> y = Variable::Make("y");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> b = MakeFloat(bstr_real, bstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	y->set_current_value<dbl>(ynum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	y->set_current_value<mpfr>(bertini::mpfr_complex(ystr_real,ystr_imag));
	
	dbl exact_dbl = xnum_dbl+anum_dbl/ynum_dbl+bnum_dbl;
	mpfr exact_mpfr{xnum_mpfr+anum_mpfr/ynum_mpfr+bnum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)+anum_dbl;
	mpfr exact_mpfr{pow(xnum_mpfr,pnum_mpfr)+anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(anum_dbl,pnum_dbl)+xnum_dbl;
	mpfr exact_mpfr{pow(anum_mpfr,pnum_mpfr)+xnum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(anum_dbl,pnum_dbl)*xnum_dbl;
	mpfr exact_mpfr{pow(anum_mpfr,pnum_mpfr)*xnum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)*anum_dbl;
	mpfr exact_mpfr{pow(xnum_mpfr,pnum_mpfr)*anum_mpfr};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl)/anum_dbl;
	mpfr exact_mpfr{pow(xnum_mpfr,pnum_mpfr)/anum_mpfr};
	
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
	bertini::mpfr_complex exact_mpfr = pow(sqrt(xnum_mpfr),anum_mpfr);

	std::shared_ptr<Variable> x = Variable::Make("x");
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = xnum_dbl/pow(anum_dbl,pnum_dbl);
	mpfr exact_mpfr{xnum_mpfr/pow(anum_mpfr,pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,anum_dbl+pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr,anum_mpfr+pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,anum_dbl*pnum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr,anum_mpfr*pnum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	std::shared_ptr<Float> p = MakeFloat(pstr_real, pstr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = pow(xnum_dbl,pnum_dbl/anum_dbl);
	mpfr exact_mpfr{pow(xnum_mpfr,pnum_mpfr/anum_mpfr)};
	
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
	mpfr exact_mpfr{sin(anum_mpfr)};
	
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
	mpfr exact_mpfr{cos(anum_mpfr)};
	
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
	mpfr exact_mpfr{tan(anum_mpfr)};
	
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
	mpfr exact_mpfr{exp(anum_mpfr)};
	
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
	mpfr exact_mpfr{sqrt(anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = sin(xnum_dbl+anum_dbl);
	mpfr exact_mpfr{sin(xnum_mpfr+anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = cos(xnum_dbl*anum_dbl);
	mpfr exact_mpfr{cos(xnum_mpfr*anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = tan(xnum_dbl/anum_dbl);
	mpfr exact_mpfr{tan(xnum_mpfr/anum_mpfr)};
	
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
	mpfr exact_mpfr{exp(-anum_mpfr)};
	
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
	
	std::shared_ptr<Variable> x = Variable::Make("x");
	std::shared_ptr<Float> a = MakeFloat(astr_real, astr_imag);
	
	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));
	
	dbl exact_dbl = sqrt(pow(xnum_dbl,anum_dbl));
	mpfr exact_mpfr{sqrt(pow(xnum_mpfr,anum_mpfr))};
	
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
	std::shared_ptr<Variable> x = Variable::Make("x");
	auto N = asin(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));

	dbl exact_dbl = asin(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr{asin(pow(xnum_mpfr,2)+bertini::mpfr_complex(1.0))};

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(arccosine_evaluate)
{
	std::shared_ptr<Variable> x = Variable::Make("x");
	auto N = acos(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));

	dbl exact_dbl = acos(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr{acos(pow(xnum_mpfr,2)+bertini::mpfr_complex(1.0))};

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(arctangent_evaluate)
{
	std::shared_ptr<Variable> x = Variable::Make("x");
	auto N = atan(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));

	dbl exact_dbl = atan(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr{atan(pow(xnum_mpfr,2)+bertini::mpfr_complex(1.0))};

	BOOST_CHECK(fabs(N->Eval<dbl>().real() - exact_dbl.real() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<dbl>().imag() - exact_dbl.imag() ) < threshold_clearance_d);
	BOOST_CHECK(fabs(N->Eval<mpfr>().real() - exact_mpfr.real() ) < threshold_clearance_mp);
	BOOST_CHECK(fabs(N->Eval<mpfr>().imag() - exact_mpfr.imag() ) < threshold_clearance_mp);

	BOOST_CHECK(!N->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(log_evaluate)
{
	std::shared_ptr<Variable> x = Variable::Make("x");
	auto N = log(pow(x,2)+1);

	x->set_current_value<dbl>(xnum_dbl);
	x->set_current_value<mpfr>(bertini::mpfr_complex(xstr_real,xstr_imag));

	dbl exact_dbl = log(pow(xnum_dbl,2)+1.0);
	mpfr exact_mpfr{log(pow(xnum_mpfr,2)+bertini::mpfr_complex(1.0))};

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
	mpfr exact_mpfr{mpfr_float("4.0")*atan(mpfr_float("1.0"))};
	
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
	mpfr exact_mpfr{exp(mpfr_float("1.0"))};
	
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
	mpfr exact_mpfr{mpfr("0.0","1.0")};
	
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
	std::shared_ptr<Variable> x = Variable::Make("x");
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

	auto rat_coeff = [](){return bertini::Rational::Make(bertini::node::Rational::Rand());};
	for (unsigned int ii=0; ii<9; ++ii)
	    polytypes[ii] = rat_coeff();

	std::shared_ptr<Node> v = rat_coeff();
	std::shared_ptr<Node> m12_2 = rat_coeff();
	std::shared_ptr<Node> cosb = rat_coeff();
	std::shared_ptr<Node> sinb = rat_coeff();
	std::shared_ptr<Node> sin2b = rat_coeff();

	std::shared_ptr<Node> N = 2*m12_2/sin2b - v*v * polytypes[9-1] / 16;

}

BOOST_AUTO_TEST_CASE(make_linear_product)
{
	using namespace bertini::node;
	
	bertini::VariableGroup vargp;
	std::shared_ptr<Variable> x = bertini::Variable::Make("x");
	std::shared_ptr<Variable> y = bertini::Variable::Make("y");
	std::shared_ptr<Variable> z = bertini::Variable::Make("z");
	std::shared_ptr<Variable> U = bertini::Variable::Make("U");
	vargp.push_back(x);
	vargp.push_back(y);
	vargp.push_back(z);
	
	// Make with automatically genearted coefficients
	std::shared_ptr<LinearProduct> linprod = LinearProduct::Make(vargp,4);
	
	// Make with user defined coefficients
	Mat<mpfr> coeff_mpfr(3,4);
	for(int ii = 0; ii < 3; ++ii)
	{
		for(int jj = 0; jj < 4; ++jj)
		{
			coeff_mpfr(ii,jj) = mpfr(ii+1, jj+1);
		}
	}
	
	std::shared_ptr<LinearProduct> linprod2 = LinearProduct::Make(vargp, coeff_mpfr);
	
}



BOOST_AUTO_TEST_CASE(eval_linear_product)
{
	using namespace bertini::node;
	
    auto x = Variable::Make("x");
    auto y = Variable::Make("y");
    auto z = Variable::Make("z");
    auto w = Variable::Make("w");
    
    
    
    VariableGroup v0{x,z,y};
    VariableGroup v1{w};
    Mat<dbl> coeff_dbl(3,4);
    Mat<mpfr> coeff_mpfr(3,4);
    
    for(int ii = 0; ii < 3; ++ii)
    {
        for(int jj = 0; jj < 4; ++jj)
        {
            coeff_dbl(ii,jj) = dbl(ii+1, jj+1);
            coeff_mpfr(ii,jj) = mpfr(ii+1, jj+1);
        }
    }
    
    std::shared_ptr<bertini::node::Node> linprod1 = bertini::LinearProduct::Make(v0, coeff_mpfr);
    
    coeff_dbl = Mat<dbl>(1,2);
    coeff_mpfr = Mat<mpfr>(1,2);
    
    for(int ii = 0; ii < 1; ++ii)
    {
        for(int jj = 0; jj < 2; ++jj)
        {
            coeff_dbl(ii,jj) = dbl(ii+3, jj+3);
            coeff_mpfr(ii,jj) = mpfr(ii+3, jj+3);
        }
    }
    
    std::shared_ptr<bertini::node::Node> linprod2 = bertini::LinearProduct::Make(v1, coeff_mpfr);
    
    
    
    std::shared_ptr<bertini::node::Node> linprod_node = (mpfr(1,1)*x + mpfr(1,2)*z + mpfr(1,3)*y+ mpfr(1,4)) * (mpfr(2,1)*x + mpfr(2,2)*z + mpfr(2,3)*y+ mpfr(2,4))*(mpfr(3,1)*x + mpfr(3,2)*z + mpfr(3,3)*y+ mpfr(3,4))*
    (mpfr(3,3)*w + mpfr(3,4));
    std::shared_ptr<bertini::node::Node> linprod = linprod1*linprod2;
    
    dbl xval_d = dbl(.5,1);
    dbl yval_d = dbl(.6,1);
    dbl zval_d = dbl(.7,1);
    dbl wval_d = dbl(2.1, -.03);
    mpfr xval_mp = mpfr(".5", "1");
    mpfr yval_mp = mpfr(".6", "1");
    mpfr zval_mp = mpfr(".7", "1");
    mpfr wval_mp = mpfr(".22", "1.11");
    
    v0[0]->set_current_value(xval_d);
    v0[1]->set_current_value(zval_d);
    v0[2]->set_current_value(yval_d);
    v0[0]->set_current_value(xval_mp);
    v0[1]->set_current_value(zval_mp);
    v0[2]->set_current_value(yval_mp);
    
    v1[0]->set_current_value(wval_d);
    v1[0]->set_current_value(wval_mp);
    
    
    dbl eval_d = linprod->Eval<dbl>();
    dbl exact_d = linprod_node->Eval<dbl>();
    mpfr eval_mp = linprod->Eval<mpfr>();
    mpfr exact_mp = linprod_node->Eval<mpfr>();
	
	
	BOOST_CHECK(fabs(eval_d.real()/exact_d.real() - 1) < threshold_clearance_d);
	BOOST_CHECK(fabs(eval_d.imag()/exact_d.imag() - 1) < threshold_clearance_d);
	BOOST_CHECK(fabs(eval_mp.real()/exact_mp.real() - 1) < threshold_clearance_mp);
	BOOST_CHECK(fabs(eval_mp.imag()/exact_mp.imag() - 1) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(linear_product_degree)
{
	using namespace bertini::node;
	
	bertini::VariableGroup vargp, vargp2, vargp3;
	std::shared_ptr<Variable> x = bertini::Variable::Make("x");
	std::shared_ptr<Variable> y = bertini::Variable::Make("y");
	std::shared_ptr<Variable> z = bertini::Variable::Make("z");
	std::shared_ptr<Variable> w = bertini::Variable::Make("w");
	vargp.push_back(x);
	vargp.push_back(y);
	vargp.push_back(z);
	
	vargp2.push_back(x);
	vargp2.push_back(z);
	vargp2.push_back(w);
	
	vargp3.push_back(w);
	
	
	
	// Make with automatically genearted coefficients
	std::shared_ptr<LinearProduct> linprod = LinearProduct::Make(vargp,4);
	
	BOOST_CHECK_EQUAL(linprod->Degree(x), 4);
	BOOST_CHECK_EQUAL(linprod->Degree(y), 4);
	BOOST_CHECK_EQUAL(linprod->Degree(z), 4);
	BOOST_CHECK_EQUAL(linprod->Degree(w), 0);
	
	BOOST_CHECK_EQUAL(linprod->Degree(vargp), 4);
	BOOST_CHECK_EQUAL(linprod->Degree(vargp2), 4);
	BOOST_CHECK_EQUAL(linprod->Degree(vargp3), 0);
	
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp)[0], 4);
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp)[1], 4);
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp)[2], 4);
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp2)[0], 4);
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp2)[2], 0);
	BOOST_CHECK_EQUAL(linprod->MultiDegree(vargp3)[0], 0);
	
}




BOOST_AUTO_TEST_CASE(linear_prod_get_linears)
{
	using namespace bertini::node;
	
    auto x = Variable::Make("x");
    auto y = Variable::Make("y");
    auto z = Variable::Make("z");
    auto w = Variable::Make("w");
    
    int num_vars = 3;
    
    
    
    VariableGroup v0{x,z,y};
    VariableGroup v1{w};
    Mat<dbl> coeff_dbl(3,4);
    Mat<mpfr> coeff_mpfr(3,4);
    
    for(int ii = 0; ii < 3; ++ii)
    {
        for(int jj = 0; jj < 4; ++jj)
        {
            coeff_dbl(ii,jj) = dbl(ii+1, jj+1);
            coeff_mpfr(ii,jj) = mpfr(ii+1, jj+1);
        }
    }
    
    std::shared_ptr<bertini::node::LinearProduct> linprod = bertini::LinearProduct::Make(v0, coeff_mpfr);
    
    
    
    std::shared_ptr<bertini::node::Node> linprod_node1 = (mpfr(1,1)*x + mpfr(1,2)*z + mpfr(1,3)*y+ mpfr(1,4));
    std::shared_ptr<bertini::node::Node> linprod_node2 = (mpfr(2,1)*x + mpfr(2,2)*z + mpfr(2,3)*y+ mpfr(2,4));
    std::shared_ptr<bertini::node::Node> linprod_node3 = (mpfr(3,1)*x + mpfr(3,2)*z + mpfr(3,3)*y+ mpfr(3,4));
    std::shared_ptr<bertini::node::Node> linprod_node12 = linprod_node1 * linprod_node3;
    std::shared_ptr<bertini::node::Node> linprod_node = linprod_node1 * linprod_node2 * linprod_node3;

    
    dbl xval_d = dbl(.5,1);
    dbl yval_d = dbl(.6,1);
    dbl zval_d = dbl(.7,1);
    dbl wval_d = dbl(2.1, -.03);
    mpfr xval_mp = mpfr(".5", "1");
    mpfr yval_mp = mpfr(".6", "1");
    mpfr zval_mp = mpfr(".7", "1");
    mpfr wval_mp = mpfr(".22", "1.11");
    
    v0[0]->set_current_value(xval_d);
    v0[1]->set_current_value(zval_d);
    v0[2]->set_current_value(yval_d);
    v0[0]->set_current_value(xval_mp);
    v0[1]->set_current_value(zval_mp);
    v0[2]->set_current_value(yval_mp);
    
    v1[0]->set_current_value(wval_d);
    v1[0]->set_current_value(wval_mp);
    
	bertini::Vec<dbl> coeff_d = linprod->GetCoeffs<dbl>(1);
	bertini::Vec<mpfr> coeff_mp = linprod->GetCoeffs<mpfr>(2);
	
	
	
	for(int ii = 0; ii < num_vars+1; ++ii)
	{
		BOOST_CHECK(fabs(coeff_d[ii].real()/coeff_dbl(1,ii).real() - 1) < threshold_clearance_d);
		BOOST_CHECK(fabs(coeff_d[ii].imag()/coeff_dbl(1,ii).imag() - 1) < threshold_clearance_d);
		BOOST_CHECK(fabs(coeff_mp[ii].real()/coeff_mpfr(2,ii).real() - 1) < threshold_clearance_mp);
		BOOST_CHECK(fabs(coeff_mp[ii].imag()/coeff_mpfr(2,ii).imag() - 1) < threshold_clearance_mp);
	}
    
    std::vector<size_t> ind{0,2};
    
//    std::shared_ptr<bertini::node::Node> linprod12 = linprod->GetLinears(ind);
//	std::shared_ptr<bertini::node::Node> linprod2 = linprod->GetLinears(1);
//    std::shared_ptr<bertini::node::Node> linprod3 = linprod->GetLinears(2);
//
//    
//    dbl evalx_d = linprod12->Eval<dbl>();
//    dbl exactx_d = linprod_node12->Eval<dbl>();
//    mpfr evalx_mp = linprod12->Eval<mpfr>();
//    mpfr exactx_mp = linprod_node12->Eval<mpfr>();
//    
//    BOOST_CHECK(fabs(evalx_d.real()/exactx_d.real() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_d.imag()/exactx_d.imag() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_mp.real()/exactx_mp.real() - 1) < threshold_clearance_mp);
//    BOOST_CHECK(fabs(evalx_mp.imag()/exactx_mp.imag() - 1) < threshold_clearance_mp);
//
//     evalx_d = linprod2->Eval<dbl>();
//     exactx_d = linprod_node2->Eval<dbl>();
//     evalx_mp = linprod2->Eval<mpfr>();
//     exactx_mp = linprod_node2->Eval<mpfr>();
//    
//    BOOST_CHECK(fabs(evalx_d.real()/exactx_d.real() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_d.imag()/exactx_d.imag() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_mp.real()/exactx_mp.real() - 1) < threshold_clearance_mp);
//    BOOST_CHECK(fabs(evalx_mp.imag()/exactx_mp.imag() - 1) < threshold_clearance_mp);
//
//     evalx_d = linprod3->Eval<dbl>();
//     exactx_d = linprod_node3->Eval<dbl>();
//     evalx_mp = linprod3->Eval<mpfr>();
//     exactx_mp = linprod_node3->Eval<mpfr>();
//    
//    BOOST_CHECK(fabs(evalx_d.real()/exactx_d.real() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_d.imag()/exactx_d.imag() - 1) < threshold_clearance_d);
//    BOOST_CHECK(fabs(evalx_mp.real()/exactx_mp.real() - 1) < threshold_clearance_mp);
//    BOOST_CHECK(fabs(evalx_mp.imag()/exactx_mp.imag() - 1) < threshold_clearance_mp);

}


BOOST_AUTO_TEST_SUITE_END()






