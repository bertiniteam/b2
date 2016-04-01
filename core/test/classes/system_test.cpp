//This file is part of Bertini 2.0.
//
//system_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  system_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

/**
\file system_test.cpp Unit testing for the bertini::System class.
*/

#include <boost/test/unit_test.hpp>



#include "bertini2/system.hpp"
#include "bertini2/system_parsing.hpp"

using System = bertini::System;
using Var = std::shared_ptr<bertini::Variable>;
using VariableGroup = bertini::VariableGroup;

using mpfr_float = bertini::mpfr_float;
using dbl = bertini::dbl;
using mpfr = bertini::mpfr;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned CLASS_TEST_MPFR_DEFAULT_DIGITS;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

BOOST_AUTO_TEST_SUITE(system_class)


/**
\class bertini::System
\test \b system_make_a_system_at_all Confirms that can default construct a System.
*/
BOOST_AUTO_TEST_CASE(system_make_a_system_at_all)
{
	System S;
}


/**
\class bertini::System
\test \b system_create_parser Confirms that can parse a system defined by a Bertini Classic style string.
*/
BOOST_AUTO_TEST_CASE(system_create_parser)
{
	System sys;
	std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";

	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();

	bertini::SystemParser<std::string::const_iterator> S;

	bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);

	BOOST_CHECK(s && iter==end);
	BOOST_CHECK(!sys.IsHomogeneous());

}





/**
\class bertini::System
\test \b system_parse_xyz_f1f2_t_pq Confirms that can parse a system defined by a Bertini Classic style string, containing two functions, powers, parameters, a path variable, and three variables.
*/
BOOST_AUTO_TEST_CASE(system_parse_xyz_f1f2_t_pq)
{

	std::string str = "variable_group x, y, z;\n function f1, f2;\n pathvariable t;\n parameter p, q;\n p = t;\n q = 1-t;\n f1 = x*y*z;\n\nf2 = p*q*x - 2^(-5);\n";




	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);

}


/**
\class bertini::System
\test \b system_parse_with_subfunctions Confirms that can parse a system with a subfunction.
*/
BOOST_AUTO_TEST_CASE(system_parse_with_subfunctions)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);

	BOOST_CHECK(!sys.IsHomogeneous());
}


/**
\class bertini::System
\test \b system_parse_around_the_unit_circle Tests parsing a system with cosine, sine, \f$\pi\f$, \f$i\f$, and subfunctions.
*/
BOOST_AUTO_TEST_CASE(system_parse_around_the_unit_circle)
{
	std::string str =
 "variable z;\nfunction H;\nparameter q1,q2;\npathvariable t;\nq1 = cos(2*Pi*(1-t));\nq2 = sin(2*Pi*(1-t));\ns = q1 + I*q2;\nH = z^2 - s;\n";



	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);

	BOOST_CHECK(!sys.IsHomogeneous());
}



/**
\class bertini::System
\test \b system_parse_around_the_unit_circle_alt Tests parsing an exponential function with time in it as a parameter.
*/
BOOST_AUTO_TEST_CASE(system_parse_around_the_unit_circle_alt)
{
	std::string str = " variable z; function H; parameter s; pathvariable t; s = exp(2*Pi*I*(1-t)); H = z^2 - s; ";


	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);

}



//TODO: uncomment this test once error handling has been done for the parsers.
// BOOST_AUTO_TEST_CASE(system_parse_x_y_not_xy)
// {
// 	std::string str = " variable x, y; function f; f = xy;";


// 	bertini::System sys;
// 	std::string::const_iterator iter = str.begin();
// 	std::string::const_iterator end = str.end();
// 	bertini::SystemParser<std::string::const_iterator> S;
// 	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
// 	BOOST_CHECK(!s && iter!=end);

// }



/**
\class bertini::System
\test \b system_differentiate_x Create a system with two functions, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	auto f1 = pow(x,2);
	auto f2 = x-1;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<dbl> v(1);
	v << 1.0;

	auto J = S.Jacobian(v);

	BOOST_CHECK_EQUAL(J(0),2.0);
	BOOST_CHECK_EQUAL(J(1),1.0);
}


/**
\class bertini::System
\test \b system_differentiate_x_and_y Create a system with two functions and two variables, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x_and_y)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> y = std::make_shared<bertini::Variable>("y");
	auto f1 = pow(x,2)*y/2;
	auto f2 = x-y;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddUngroupedVariable(y);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<dbl> v(2);
	v << 1.0 , 2.0;

	auto J = S.Jacobian(v);

	BOOST_CHECK_THROW(S.Jacobian(v,dbl(0.5)), std::runtime_error);
}


/**
\class bertini::System
\test \b system_differentiate_x_and_t Create a system with two functions, one variable, one time variable, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x_and_t)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");
	auto f1 = (1-t)*x + t*(1-x);
	auto f2 = x-t;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<dbl> v(1);
	v << 1.0;
	dbl time(0.5,0.2);
	auto J = S.Jacobian(v,time);

	BOOST_CHECK_THROW(S.Jacobian(v), std::runtime_error);
}




/**
\class bertini::System
\test \b system_homogenize_multiple_variable_groups Homogenize a system with multiple variable groups.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_multiple_variable_groups)
{
	Var x1 = std::make_shared<bertini::Variable>("x1");
	Var x2 = std::make_shared<bertini::Variable>("x2");

	Var y1 = std::make_shared<bertini::Variable>("y1");
	Var y2 = std::make_shared<bertini::Variable>("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*x1 + y1*y2 + x1 + y2 - 1;


	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1);
	S.AddFunction(f2);

	BOOST_CHECK(!S.IsHomogeneous());

	BOOST_CHECK_EQUAL(S.NumHomVariables(),0);

	S.Homogenize();

	BOOST_CHECK(S.IsHomogeneous());
}



/**
\class bertini::System
\test \b system_reorder_by_degree_decreasing For a system with four functions, re-order the functions so they are in decreasing degree.
*/
BOOST_AUTO_TEST_CASE(system_reorder_by_degree_decreasing)
{
	Var x1 = std::make_shared<bertini::Variable>("x1");
	Var x2 = std::make_shared<bertini::Variable>("x2");

	Var y1 = std::make_shared<bertini::Variable>("y1");
	Var y2 = std::make_shared<bertini::Variable>("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*pow(x1,2) + y1*y2 + x1 + y2 - 1;

	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1); // deg 2
	S.AddFunction(f2); // deg 3 
	S.AddFunction(pow(x1,4) + pow(y2,5)); // deg 5
	S.AddFunction(x1 + x2 + y1 + y2); // deg 1


	BOOST_CHECK(!S.IsHomogeneous());
	S.ReorderFunctionsByDegreeDecreasing();

	auto degs = S.Degrees();

	for (auto d = degs.begin(); d != degs.end()-1; d++)
	{
		BOOST_CHECK(*d >= *(d+1));
	}

}



/**
\class bertini::System
\test \b system_reorder_by_degree_increasing For a system with four functions, re-order the functions so they are in increasing degree.
*/
BOOST_AUTO_TEST_CASE(system_reorder_by_degree_increasing)
{
	Var x1 = std::make_shared<bertini::Variable>("x1");
	Var x2 = std::make_shared<bertini::Variable>("x2");

	Var y1 = std::make_shared<bertini::Variable>("y1");
	Var y2 = std::make_shared<bertini::Variable>("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*pow(x1,2) + y1*y2 + x1 + y2 - 1;


	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1);
	S.AddFunction(f2);

	BOOST_CHECK(!S.IsHomogeneous());
	S.ReorderFunctionsByDegreeIncreasing();

	auto degs = S.Degrees();

	for (auto d = degs.begin(); d != degs.end()-1; d++)
	{
		BOOST_CHECK(*d <= *(d+1));
	}

}





/**
\class bertini::System
\test \b system_evaluate_double Evaluate a system in double precision.
*/
BOOST_AUTO_TEST_CASE(system_evaluate_double)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);

	Vec<dbl> values(2);

	values(0) = dbl(2.0);
	values(1) = dbl(3.0);

	Vec<dbl> v = sys.Eval(values);

	BOOST_CHECK_EQUAL(v(0), 36.0);


	auto J = sys.Jacobian(values);

	double x1 = 2;
	double x2 = 3;

	BOOST_CHECK_EQUAL(J(0,0), 2*x1*x2*x2);
	BOOST_CHECK_EQUAL(J(0,1), x1*x1*2*x2);
}


/**
\class bertini::System
\test \b system_evaluate_mpfr Evaluate a system in multiple precision.
*/
BOOST_AUTO_TEST_CASE(system_evaluate_mpfr)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);

	Vec<mpfr> values(2);

	values(0) = mpfr(2.0);
	values(1) = mpfr(3.0);

	Vec<mpfr> v = sys.Eval(values);

	BOOST_CHECK_EQUAL(v(0), mpfr(36.0));


	auto J = sys.Jacobian(values);

	mpfr x1 = 2;
	mpfr x2 = 3;
	
	BOOST_CHECK_EQUAL(J(0,0), mpfr(2.0)*x1*x2*x2);
	BOOST_CHECK_EQUAL(J(0,1), x1*x1*mpfr(2.0)*x2);
}




/**
\class bertini::System
\test \b add_two_systems Test the arithmetic sum of two Systems.
*/
BOOST_AUTO_TEST_CASE(add_two_systems)
{
	bertini::System sys1, sys2;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys1.AddVariableGroup(vars);  
	sys1.AddFunction(y+1);
	sys1.AddFunction(x*y);

	sys2.AddVariableGroup(vars);  
	sys2.AddFunction(-y-1);
	sys2.AddFunction(-x*y);

	sys1+=sys2;


	Vec<dbl> values(2);

	values << dbl(2.0), dbl(3.0);

	auto v = sys1.Eval(values);

	BOOST_CHECK_EQUAL(v(0), 0.0);
	BOOST_CHECK_EQUAL(v(1), 0.0);

	auto deg = sys1.Degrees();

	BOOST_CHECK_EQUAL(deg.size(),2);
	if (deg.size()==2)
	{
		BOOST_CHECK_EQUAL(deg[0],1);
		BOOST_CHECK_EQUAL(deg[1],2);
	}


}



/**
\class bertini::System
\test \b system_differentiate_wrt_time_linear Test the arithmetic sum of two Systems.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_wrt_time_linear)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");
	auto f1 = (1-t)*x + t*(1-x);
	auto f2 = x-t;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<dbl> v(1);
	v << 1.0;
	dbl time(0.5,0.2);
	auto dS_dt = S.TimeDerivative(v,time);

	BOOST_CHECK( abs(dS_dt(0) - dbl(-1) ) < threshold_clearance_d);
	BOOST_CHECK( abs(dS_dt(1) - dbl(-1) ) < threshold_clearance_d);


}





/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_aff_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_aff_group)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	VariableGroup vars{x, y};
	sys.AddVariableGroup(vars);

	sys.Homogenize();

	Vec<dbl> v(3);
	v << dbl(2,3), dbl(3,4), dbl(4,5);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),2);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_dehomogenize_FIFO_two_aff_groups Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_two_aff_groups)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z"), w = std::make_shared<bertini::Variable>("w");
	VariableGroup vars{x, y};
	VariableGroup vars2{z, w};
	sys.AddVariableGroup(vars);
	sys.AddVariableGroup(vars2);

	sys.Homogenize();

	Vec<dbl> v(6);
	v << dbl(2,3), dbl(3,4), dbl(4,5), 
		 dbl(5,6), dbl(6,7), dbl(7,8);

	auto d = sys.DehomogenizePoint(v);


	BOOST_CHECK_EQUAL(d.size(),4);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(4)/v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(5)/v(3)) < threshold_clearance_d);
}



/**
\class bertini::System
\test \b system_dehomogenize_FIFO_two_aff_groups_one_hom_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_two_aff_groups_one_hom_group)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z"), w = std::make_shared<bertini::Variable>("w");
	Var h1 = std::make_shared<bertini::Variable>("h1"), h2 = std::make_shared<bertini::Variable>("h2");
	VariableGroup vars{x, y};
	VariableGroup vars2{h1,h2};
	VariableGroup vars3{z, w};
	sys.AddVariableGroup(vars);
	sys.AddHomVariableGroup(vars2);
	sys.AddVariableGroup(vars3);

	sys.Homogenize();

	Vec<dbl> v(8);
	v << dbl(2,3), dbl(3,4), dbl(4,5), 
		 dbl(10,11), dbl(11,12),
		 dbl(5,6), dbl(6,7), dbl(7,8);

	auto d = sys.DehomogenizePoint(v);


	BOOST_CHECK_EQUAL(d.size(),6);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(4)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(4) - v(6)/v(5)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(5) - v(7)/v(5)) < threshold_clearance_d);
}



/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_hom_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_hom_group)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	VariableGroup vars{x, y};
	sys.AddHomVariableGroup(vars);

	sys.Homogenize();

	Vec<dbl> v(2);
	v << dbl(2,3), dbl(3,4);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),2);

	BOOST_CHECK(abs(d(0) - v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(1)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_hom_group_two_ungrouped_vars Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_hom_group_two_ungrouped_vars)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z"), w = std::make_shared<bertini::Variable>("w");
	VariableGroup vars{x, y};

	sys.AddHomVariableGroup(vars);
	sys.AddUngroupedVariable(z);
	sys.AddUngroupedVariable(w);

	sys.Homogenize();

	Vec<dbl> v(4);
	v << dbl(2,3), dbl(3,4), dbl(4,5), dbl(5,6);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),4);

	BOOST_CHECK(abs(d(0) - v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(2) - v(2)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(3)) < threshold_clearance_d);
}

/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_aff_group_two_ungrouped_vars_another_aff_grp_hom_grp Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_aff_group_two_ungrouped_vars_another_aff_grp_hom_grp)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z"), w = std::make_shared<bertini::Variable>("w");
	Var h1 = std::make_shared<bertini::Variable>("h1"), h2 = std::make_shared<bertini::Variable>("h2");
	Var u1 = std::make_shared<bertini::Variable>("u1"), u2 = std::make_shared<bertini::Variable>("u2");

	VariableGroup vars{x,y};
	VariableGroup vars2{z,w};

	VariableGroup vars3{h1,h2};

	sys.AddVariableGroup(vars);
	sys.AddUngroupedVariable(u1);
	sys.AddUngroupedVariable(u2);
	sys.AddVariableGroup(vars2);
	sys.AddHomVariableGroup(vars3);

	sys.Homogenize();

	Vec<dbl> v(10);
	v << dbl(2,3), dbl(3,4), dbl(4,5), 
		 dbl(10,11), dbl(11,12),
		 dbl(5,6), dbl(6,7), dbl(7,8),
		 dbl(12,13), dbl(13,14);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),8);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(4)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(4) - v(6)/v(5)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(5) - v(7)/v(5)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(6) - v(8)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(7) - v(9)) < threshold_clearance_d);
}










/**
\class bertini::System
\test \b system_estimate_coeff_bound_linear Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_linear)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction((1-t)*x + t*(1-x));
	S.AddFunction(x-t);

	mpfr_float coefficient_bound = S.CoefficientBound();
	BOOST_CHECK(coefficient_bound < mpfr_float("3"));
	BOOST_CHECK(coefficient_bound > mpfr_float("0.5"));
}


/**
\class bertini::System
\test \b system_estimate_coeff_bound_quartic Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_quartic)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	mpfr_float coefficient_bound = sys.CoefficientBound();
	BOOST_CHECK(coefficient_bound < mpfr_float("5"));
	BOOST_CHECK(coefficient_bound > mpfr_float("2"));
}



/**
\class bertini::System
\test \b system_estimate_coeff_bound_quartic Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_homogenized_quartic)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	sys.Homogenize();
	sys.AutoPatch();

	mpfr_float coefficient_bound = sys.CoefficientBound();
	BOOST_CHECK(coefficient_bound < mpfr_float("10"));
	BOOST_CHECK(coefficient_bound > mpfr_float("2"));
}

/**
\class bertini::System
\test \b system_estimate_degree_bound_linear Test the estimation of the degree in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_degree_bound_linear)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction((1-t)*x + t*(1-x));
	S.AddFunction(x-t);

	mpfr_float degree_bound = S.DegreeBound();
	BOOST_CHECK(degree_bound == mpfr_float("1"));
}


/**
\class bertini::System
\test \b system_estimate_degree_bound_quartic Test the estimation of the degree in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_degree_bound_quartic)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	mpfr_float degree_bound = sys.DegreeBound();
	BOOST_CHECK(degree_bound == mpfr_float("4"));
}




/**
\class bertini::System
\test \b system_multiply_by_node Ensure that multiplication of a system by a node doesn't affect other copies of a system
*/
BOOST_AUTO_TEST_CASE(system_multiply_by_node)
{
	bertini::System sys1, sys2;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars{x,y,z};

	sys1.AddVariableGroup(vars);  
	sys1.AddFunction(x);
	sys1.AddFunction(y);
	sys1.AddFunction(z);

	sys2.AddVariableGroup(vars);  
	sys2.AddFunction(y+x*y + mpfr_float("0.5"));
	sys2.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys2.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	Var t = std::make_shared<bertini::node::Variable>("t");

	auto sys_copy1 = t*sys1;
	auto sys_copy2 = (1-t)*sys2;

	auto sys_copy3 = sys_copy1 + sys_copy2; // this line couples the two systems...  this coupling is total garbage!!!  this 'homotopy' should never be used.
}








BOOST_AUTO_TEST_SUITE_END()




