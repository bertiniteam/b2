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



#include <boost/test/unit_test.hpp>



#include "system.hpp"
#include "system_parsing.hpp"

using System = bertini::System;
template<typename T>
	using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	template<typename T>
	using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

BOOST_AUTO_TEST_SUITE(system_class)


BOOST_AUTO_TEST_CASE(system_make_a_system_at_all)
{
	System S;
}

BOOST_AUTO_TEST_CASE(system_create_parser)
{
	System sys;
	std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";

	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();


	bertini::SystemParser<std::string::const_iterator> S;


	bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);

	std::cout << std::string(iter,end);


	BOOST_CHECK(s && iter==end);
}






BOOST_AUTO_TEST_CASE(system_parse_xyz_f1f2_t_pq)
{

	std::string str = "variable_group x, y, z;\n function f1, f2;\n pathvariable t;\n parameter p, q;\n p = t;\n q = 1-t;\n f1 = x*y*z;\n\nf2 = p*q*x - 2^(-5);\n";




	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);


//	std::cout << str << " got parsed into the following system in memory:\n\n";
//
//	std::cout << sys << std::endl;
//
//
//	bertini::Vec<double> var_values(3);
//	var_values << 1, 2, 3;
//
//	double t = 2.0/3.0;
//	bertini::Vec<double> sys_values(2);
//
//	sys_values = sys.Eval(var_values,t);
//
//	std::cout << "system evaluated at\n" << var_values << " is " << sys_values << "\n";

}



BOOST_AUTO_TEST_CASE(system_parse_with_subfunctions)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(s && iter==end);
}


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
}




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




BOOST_AUTO_TEST_CASE(system_parse_x_y_not_xy)
{
	std::string str = " variable x, y; function f; f = xy;";


	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	BOOST_CHECK(!s && iter!=end);

}




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

	std::cout << S << std::endl;
	std::cout << "v\n" << v << std::endl;
	std::cout << "J\n" << J << std::endl;
}



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

	std::cout << S << std::endl;
	std::cout << "v\n" << v << std::endl;
	std::cout << "J\n" << J << std::endl;
}




BOOST_AUTO_TEST_SUITE_END()




