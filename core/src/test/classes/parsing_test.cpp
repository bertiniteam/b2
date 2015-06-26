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

//  parsing_test.cpp
//
//  Created by Collins, James B. on 6/9/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "bertini.hpp"
#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>



#include <boost/test/unit_test.hpp>

using IteratorType = std::string::const_iterator;

using VariableParser = bertini::VariableParser<IteratorType>;
using FunctionParser = bertini::FunctionParser<IteratorType>;
using Variable = bertini::Variable;
using boost::spirit::ascii::space;
using bertini::qi;

using mpfr_float = boost::multiprecision::mpfr_float;



BOOST_AUTO_TEST_SUITE(parsing_class)

/////////// Parsing Variables(Should Work) ///////////////////

BOOST_AUTO_TEST_CASE(construct_variable_parser){
    VariableParser variable_parser;
    
    
//    BOOST_CHECK(abs(N->Eval<mpfr>().imag() - mpfr_float("38.08")) < 1e-15);
    
}


BOOST_AUTO_TEST_CASE(parse_variable_x){
    VariableParser variable_parser;
    
    std::string str = "x";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}


BOOST_AUTO_TEST_CASE(parse_variable_named_myvariable){
    VariableParser variable_parser;
    
    std::string str = "myvariable";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}



BOOST_AUTO_TEST_CASE(parse_variable_x_y_z){
    VariableParser variable_parser;
    
    std::string str = "x,y,  z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}


BOOST_AUTO_TEST_CASE(parse_variable_x_y_z_with_numbers){
    VariableParser variable_parser;
    
    std::string str = "x2, y3,  z56";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}



BOOST_AUTO_TEST_CASE(parse_variable_three_with_long_names){
    VariableParser variable_parser;
    
    std::string str = "dog2, giraffe45, kettle98747";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}


BOOST_AUTO_TEST_CASE(parse_variable_maple_x_0_1_2){
    VariableParser variable_parser;
    
    std::string str = "x[0], x[108], x[12]";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(ret && iter == end);
}


/////////// Parsing Variables(Should NOT Work) ///////////////////
BOOST_AUTO_TEST_CASE(parse_variable_maple_x_no_number){
    VariableParser variable_parser;
    
    std::string str = "x[]";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
   
    
    BOOST_CHECK(!(ret && iter == end));
}


BOOST_AUTO_TEST_CASE(parse_variable_number_first){
    VariableParser variable_parser;
    
    std::string str = "y, 2x";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(!(ret && iter == end));
}


BOOST_AUTO_TEST_CASE(parse_variable_underscore_first){
    VariableParser variable_parser;
    
    std::string str = "y, _x";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    BOOST_CHECK(!(ret && iter == end));
}












////////////Parsing Functions (Should Work)//////////////
BOOST_AUTO_TEST_CASE(parse_pi){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    FunctionParser func_parser(&junk, variable_parser.variable());
    
    str = "pi";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (3.14159265358979323846264338328)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("3.14159265358979323846264338328")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("0")) < 1e-15);
    
    str = "Pi";
    iter = str.begin();
    end = str.end();
    
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    N.SetRoot(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (3.14159265358979323846264338328)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("3.14159265358979323846264338328")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("0")) < 1e-15);

}


BOOST_AUTO_TEST_CASE(parse_e){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    FunctionParser func_parser(&junk, variable_parser.variable());
    
    str = "e";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (2.71828182845904523536028747135)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("2.71828182845904523536028747135")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("0")) < 1e-15);
    
    str = "E";
    iter = str.begin();
    end = str.end();
    
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    N.SetRoot(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (2.71828182845904523536028747135)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("2.71828182845904523536028747135")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("0")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_i){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    FunctionParser func_parser(&junk, variable_parser.variable());
    
    str = "i";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (0)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (1.0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("0.0")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("1.0")) < 1e-15);
    
    str = "I";
    iter = str.begin();
    end = str.end();
    
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    N.SetRoot(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (0)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (1.0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("0.0")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("1.0")) < 1e-15);
    
    str = "1i";
    iter = str.begin();
    end = str.end();
    
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    N.SetRoot(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (0)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (1.0)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("0.0")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("1.0")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_number){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    FunctionParser func_parser(&junk, variable_parser.variable());
    
    str = "3.4+i*5.6";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (3.4)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (5.6)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("3.4")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("5.6")) < 1e-15);
}



BOOST_AUTO_TEST_CASE(parse_func_x_plus_num_plus_y_plus_z){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    junk[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    junk[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    junk[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    junk[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    junk[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    junk[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));

    FunctionParser func_parser(&junk, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "x+"+numstr + "+y+z";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (12.9)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (12.1)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("12.9")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("12.1")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_x_times_y_times_num_times_z){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    
    
    std::vector<std::shared_ptr<Variable> > junk;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        junk.push_back(std::make_shared<Variable>());
    }
    junk[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    junk[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    junk[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    junk[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    junk[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    junk[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&junk, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "x*y*"+numstr + "*z";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (2149.6134)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (2777.7882)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("2149.6134")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("2777.7882")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_linear_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "5.3*x+8.5*y-i*4*z+"+numstr;
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (64.63)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (121.08)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("64.63")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("121.08")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_degree_2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "x^2+y^2+z^2+x*y+x*z+y*z+x+y+z+"+numstr;
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-2.23)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (191.68)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-2.23")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("191.68")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_prod_of_linears){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "(x-2)*(y-3)*(z+9)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (34.239)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (485.397)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("34.239")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("485.397")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_prod_of_linears_with_powers){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "(x-2)^3*(y-3)^2*(z+9)^5";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-1.904753713e8)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (9.840123228e8)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-1.904753713631640354e8")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("9.840123227418186978e8")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_degree2_poly_squared){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "(x^2+4*y^2-7*x*y+5.4*z-6.8)^2";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-45747.7152)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-6893.7086)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-45747.7152")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-6893.7086")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_degree2_poly_over_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "(x^2+4*y^2-7*x*y+5.4*z-6.8)/(z^2-4*y*z+x-1i*4+x^2)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (0.380706454830044499787591749666)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-.203586886888085377444986761112)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("0.380706454830044499787591749666")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-.203586886888085377444986761112")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_negate_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,4.1));
    variables[1]->set_current_value<dbl>(std::complex<double>(8.8,9.9));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,-7.5));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","4.1"));
    variables[1]->set_current_value<mpfr>(bertini::complex("8.8","9.9"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","-7.5"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "-(x^2+4*y^2-7*x*y+5.4*z-6.8)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (16.07)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-214.49)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("16.07")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-214.49")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_sine_of_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,0));
    variables[1]->set_current_value<dbl>(std::complex<double>(0,3.7));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,0));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","0"));
    variables[1]->set_current_value<mpfr>(bertini::complex("0","3.7"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","0"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "sin(x^2-y*z+1+z^2-6*y+x*y)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-2.00309747340660133449118482520)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (2.44611810732349423982109591283)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-2.00309747340660133449118482520")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("2.44611810732349423982109591283")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_cosine_of_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,0));
    variables[1]->set_current_value<dbl>(std::complex<double>(0,3.7));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,0));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","0"));
    variables[1]->set_current_value<mpfr>(bertini::complex("0","3.7"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","0"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "cos(x^2-y*z+1+z^2-6*y+x*y)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-2.57013763905097552259551835195)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-1.90643992212148106701042535573)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-2.57013763905097552259551835195")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-1.90643992212148106701042535573")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_tangent_of_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,0));
    variables[1]->set_current_value<dbl>(std::complex<double>(0,3.7));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,0));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","0"));
    variables[1]->set_current_value<mpfr>(bertini::complex("0","3.7"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","0"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "tan(x^2-y*z+1+z^2-6*y+x*y)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (0.0473489535082370108790555686969)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-.986867786387826653453823677439)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("0.0473489535082370108790555686969")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-.986867786387826653453823677439")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_exp_of_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,0));
    variables[1]->set_current_value<dbl>(std::complex<double>(0,3.7));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,0));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","0"));
    variables[1]->set_current_value<mpfr>(bertini::complex("0","3.7"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","0"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "exp(x^2-y*z+1+z^2-6*y+x*y)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    std::cout << std::string(iter,end) << std::endl;
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (-3.54539405822929293604828762125e6)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-1.23665457403576246601391981124e7)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("-3.54539405822929293604828762125e6")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-1.23665457403576246601391981124e7")) < 1e-15);
}


BOOST_AUTO_TEST_CASE(parse_func_sqrt_of_degree2_poly){
    VariableParser variable_parser;
    
    std::string str = "x,y,z";
    IteratorType iter = str.begin();
    IteratorType end = str.end();
    
    bool ret = qi::phrase_parse(iter,end,variable_parser.start(),space);
    
    
    std::vector<std::shared_ptr<Variable> > variables;
    for(int ii = 0; ii < variable_parser.var_count(); ++ii)
    {
        variables.push_back(std::make_shared<Variable>());
    }
    variables[0]->set_current_value<dbl>(std::complex<double>(3.1,0));
    variables[1]->set_current_value<dbl>(std::complex<double>(0,3.7));
    variables[2]->set_current_value<dbl>(std::complex<double>(-2.4,0));
    variables[0]->set_current_value<mpfr>(bertini::complex("3.1","0"));
    variables[1]->set_current_value<mpfr>(bertini::complex("0","3.7"));
    variables[2]->set_current_value<mpfr>(bertini::complex("-2.4","0"));
    
    FunctionParser func_parser(&variables, variable_parser.variable());
    
    std::string numstr = "(3.4+i*5.6)";
    str = "sqrt(x^2-y*z+1+z^2-6*y+x*y)";
    iter = str.begin();
    end = str.end();
    
    std::shared_ptr<Node> node;
    ret = qi::phrase_parse(iter,end,func_parser.root_rule_,space, node);
    bertini::Function N(node);
    
    BOOST_CHECK(ret && iter == end);
    BOOST_CHECK(fabs(N.Eval<dbl>().real() - (4.05241928222456925923496781407)) < 1e-15);
    BOOST_CHECK(fabs(N.Eval<dbl>().imag() - (-.228258710557764074535887077743)) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().real() - mpfr_float("4.05241928222456925923496781407")) < 1e-15);
    BOOST_CHECK(abs(N.Eval<mpfr>().imag() - mpfr_float("-.228258710557764074535887077743")) < 1e-15);
}






BOOST_AUTO_TEST_SUITE_END()
