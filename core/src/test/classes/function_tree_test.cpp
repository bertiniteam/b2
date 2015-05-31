//
//  function_tree_test.cpp
//  b2Test
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <vector>


#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>

namespace qi = boost::spirit::qi;

////////// TESTING /////////////
int ManualConstruction();
int InteractiveTest();
int FixedInputTest();
int FixedInputTest2();


int Node::tabcount = 0;


// set the values of the variables, by asking the user for their values.
void set_variable_values(std::vector<std::shared_ptr<Variable> > & variables);

// print variable values to the screen
void print_variable_values(const std::vector<std::shared_ptr<Variable> > & variables);

// get the names of the variables from the user.
std::vector<std::shared_ptr<Variable> > get_var_group_from_user(qi::symbols<char,int> *variable_names);


std::shared_ptr<Node> get_function_from_user(const std::vector<std::shared_ptr<Variable> > & variables, qi::symbols<char,int>  variable_names);

////////// TESTING /////////////




int main()
{
//	FixedInputTest();
//	FixedInputTest2();
	ManualConstruction();
//	InteractiveTest();

	
	return 0;
}




int ManualConstruction(){
	
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
	std::shared_ptr<Variable> y = std::make_shared<Variable>("y");
	
	std::shared_ptr<Node> N = pow(pow(x+x*y,2),2);
//
//	std::shared_ptr<Node> N = x;
	N *= N;

	
	x->set_current_value(3.0);
	y->set_current_value(4.0);
	
	
	std::cout << "N = " << *N << "\n";
	N->Reset();
	
	std::cout << N->Eval<dbl>() << std::endl;
	
	return 0;
}




int FixedInputTest()
{
	
	
	
	rand(); // burn one for good measure
	
	std::string variable_string = "x, y";  // parse it up, boyee

	
	std::vector<std::shared_ptr<Variable> > variables; // this is produced by the variable parser
	qi::symbols<char,int> variable_names; // these will be made by the variable parser
	
	
	bertini::VariableParser<std::string::const_iterator> var_parser; //create a parser which will be used to parse the line containing the names of the variables.
	
	std::string::const_iterator iter = variable_string.begin(); // point to where you want to start parsing
	std::string::const_iterator end = variable_string.end();    // where you want to end parsing.
	
	// actually do the parsing
	qi::phrase_parse(iter,end,var_parser.start(),boost::spirit::ascii::space);
	
	variable_names = var_parser.variable();
	variables = var_parser.get_var_group();
	
	
	
	
	
	// valid strings
//	std::string function_string = "x^2^2";
//	std::string function_string = "x*x*x + y*y*y / x^2";
//	std::string function_string = "1E10";
//	std::string function_string = "x*y";
//	std::string function_string = "x+1";
//	std::string function_string = "-y^-x";
//	std::string function_string = "1^(x*y^2 + 1)";
	
	
	// invalid strings
	std::string function_string = "z";
	
	iter = function_string.begin();
	end  = function_string.end();
	
	
	bertini::BrakeParser   <std::string::const_iterator> brakerulez(&variables, variable_names);// ;
	std::shared_ptr<Node> F;
	bool succeeded = phrase_parse(iter, end, brakerulez,boost::spirit::ascii::space, F);
	
	if (! (succeeded && iter==end) ) {
		std::cout << "-------------------------\n";
		std::cout << "Parsing failed\n";
		std::cout << "stopped at: \": " << std::string(iter, end) << "\"\n";
		std::cout << "-------------------------\n";
	}
	
	
	
	for (auto iter : variables) {
		iter->set_current_value(double(rand())/RAND_MAX);
	}
	
	
	if (!F) {
		std::cout << "F is not populated for some reason" << std::endl;
	}
	else
	{
		std::cout << "\n\nF = "  << *F << std::endl;
		std::cout << "F(";
		print_variable_values(variables);
		std::cout << ") = " << std::endl;
		std::cout << F->Eval<dbl>();
	}
	
	
	return 0;
}





int FixedInputTest2()
{
	
	
	rand(); // burn one for good measure
	
	std::string variable_string = "x_1, x_2";  // parse it up, boyee
	
	
	std::vector<std::shared_ptr<Variable> > variables; // this is produced by the variable parser
	qi::symbols<char,int> variable_names; // these will be made by the variable parser
	
	
	bertini::VariableParser<std::string::const_iterator> var_parser; //create a parser which will be used to parse the line containing the names of the variables.
	
	std::string::const_iterator iter = variable_string.begin(); // point to where you want to start parsing
	std::string::const_iterator end = variable_string.end();    // where you want to end parsing.
	
	// actually do the parsing
	auto variable_parsing_succeeded = qi::phrase_parse(iter,end, var_parser.start(), boost::spirit::ascii::space);
	
	if (!(iter==end && variable_parsing_succeeded)) {
		std::cout << "variable parsing failed." << std::endl;
		std::cout << std::string(iter,end);
		return 1;
	}
	
	variable_names = var_parser.variable();
	variables = var_parser.get_var_group();
	
	for (auto iter : variables) {
		std::cout << *iter << std::endl;
	}
	
	
	
	// valid strings
		std::string function_string = "x_1^2^x_1";
	//	std::string function_string = "x*x*x + y*y*y / x^2";
	//	std::string function_string = "1E10";
	//	std::string function_string = "x*y";
	//	std::string function_string = "x+1";
	//	std::string function_string = "-y^-x";
	//	std::string function_string = "1^(x*y^2 + 1)";
	
	
	iter = function_string.begin();
	end  = function_string.end();
	
	
	bertini::BrakeParser   <std::string::const_iterator> brakerulez(&variables, variable_names);// ;
	std::shared_ptr<Node> F;
	bool succeeded = phrase_parse(iter, end, brakerulez,boost::spirit::ascii::space, F);
	
	if (! (succeeded && iter==end) ) {
		std::cout << "-------------------------\n";
		std::cout << "Parsing failed\n";
		std::cout << "stopped at: \": " << std::string(iter, end) << "\"\n";
		std::cout << "-------------------------\n";
	}
	
	
	
	for (auto iter : variables) {
		iter->set_current_value(double(rand())/RAND_MAX);
	}
	
	
	if (!F) {
		std::cout << "F is not populated for some reason" << std::endl;
	}
	else
	{
		std::cout << "\n\nF = "  << *F << std::endl;
		std::cout << "F(";
		print_variable_values(variables);
		std::cout << ") = " << std::endl;
		std::cout << F->Eval<dbl>();
	}
	
	
	return 0;
}





int InteractiveTest()
{
	
	
	std::cout.precision(10);
	
	// 1. Read in the variables //
	std::cout << "/////////////////////////////////////////////////////////\n\n";
	std::cout << "\t\tPolynomial Playground\n\n";
	std::cout << "/////////////////////////////////////////////////////////\n\n";
	
	
	//    typedef parser::variables<iterator_type> variables;
	//    typedef parser::polynomial<iterator_type> polys;
	
	
	
	
	
	
	// 2. Main loop //
	
	
	
	
	std::vector<std::shared_ptr<Variable> > variables;
//	Node* base_node{nullptr};
	
	qi::symbols<char,int> variable_names;
	
	std::shared_ptr<Node> F;
	
	bool cont = true;
	while (cont)
	{
		
		std::cout << "/////////////////////////////////////////////////////////\n\n";
		std::cout << "[v] - declare variables (resets polynomial) ..or\n";
		std::cout << "[n] - input a new polynomial ...or\n";
		std::cout << "[e] - evaluate a polynomial...or\n";
		std::cout << "[s] - set the variable values...or\n";
		std::cout << "[c] - clear all data, start over ...or\n";
		std::cout << "[q] - quit\n\n";
		std::cout << ">> ";
		std::string str;
		std::getline(std::cin,str);
		char choice = tolower(str[0]);
		
		switch (choice) {
			case 'v':
			{
				if (F==nullptr) {
					F.reset();
				}
				variables = get_var_group_from_user(&variable_names);
				break;
			}
			case 's':
			{
				set_variable_values(variables);
				if (F!=nullptr) {
					F->Reset();
				}
				
				break;
			}
			case 'e':
			{
				if(F != nullptr){
					std::cout << "f = " << *F << std::endl;
					std::cout << "test_function(";
					print_variable_values(variables);
					std::cout << ")  =  " << F->Eval<dbl>() << std::endl;
				}
				else
					std::cout << "No function defined!\n";
				break;
			}
				
			case 'n':
			{
				if (variables.size()>0){
//					if (base_node==nullptr) {
//						delete base_node;
//						base_node = nullptr;
//					}
					F = get_function_from_user(variables,variable_names);
				}
				else
					std::cout << "no variables declared.  declare variables first." << std::endl;
				break;
			}
			
			case 'c':
			{
				F.reset();
//				if (base_node==nullptr) {
//					delete base_node;
//					base_node = nullptr;
//				}
				variables.clear();
				break;
			}
				
			case 'q':
			{
				cont = false;
				break;
			}
			default:
			{
				std::cout << "invalid entry: " << choice << std::endl;
				break;
			}
		}
		

	}
	
	return 0;
}



void set_variable_values(std::vector<std::shared_ptr<Variable> > & variables)
{
	std::cout << "enter the values of the variables, as prompted\n\n";
	for (auto iter=variables.begin(); iter!=variables.end(); iter++) {
		std::cout << (*iter)->name() << " = ";
		std::string str;
		std::getline(std::cin,str);
		(*iter)->set_current_value(std::stod(str));
	}
}




void print_variable_values(const std::vector<std::shared_ptr<Variable> > & variables)
{
	for (auto iter=variables.begin(); iter!=variables.end(); iter++) {
		std::cout << (*iter)->name() << " = " << (*iter)->Eval<double>();
		if (iter!=variables.end()-1) {
			std::cout << ", ";
		}
	}
}


std::vector<std::shared_ptr<Variable> > get_var_group_from_user(qi::symbols<char,int> *variable_names)
{
	
	
	
	bertini::VariableParser<std::string::const_iterator> var_parser; //create a parser which will be used to parse the line containing the names of the variables.
	
	
	std::cout << "Enter the names of the variables separated by commas ...\n\n";
	
	std::string str;
	std::getline(std::cin, str);
	
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	
	qi::phrase_parse(iter,end,var_parser.start(),boost::spirit::ascii::space);
	
	*variable_names = var_parser.variable();
	
	// get the vector from the parser.
	return var_parser.get_var_group();
}





std::shared_ptr<Node> get_function_from_user(const std::vector<std::shared_ptr<Variable> > & variables, qi::symbols<char,int> variable_names)
{
	
	
	std::cout << "enter your polynomial.  the variables are:\n\t\t";
	for (auto iter=variables.begin(); iter!=variables.end(); iter++) {
		std::cout << (*iter)->name();
		if (iter!=variables.end()-1) {
			std::cout << ", ";
		}
	}
	std::cout << "\n\n>> ";
	
	std::string str;
	std::getline(std::cin, str);
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	
	
	bertini::BrakeParser   <std::string::const_iterator> brakerulez(&variables, variable_names);// ;
	std::shared_ptr<Node> F;
	bool s = phrase_parse(iter, end, brakerulez,boost::spirit::ascii::space, F);
	
	if (! (s && iter==end) ) {
		std::cout << "parser failed.  " << std::endl;
		std::string rest(iter, end);
		std::cout << "-------------------------\n";
		std::cout << "Parsing failed\n";
		std::cout << "stopped at: \": " << rest << "\"\n";
		std::cout << "-------------------------\n";
	}
	
	return F;
	
//	FunctionParser<std::string::const_iterator> func_parser(&variables, variable_names);
//	Node *base_node{nullptr};
//	//        test_function.entry_node().reset();
//	bool r = phrase_parse(iter, end, func_parser,boost::spirit::ascii::space, base_node);
//	//        test_function.AddChild(std::move(parse_ret) );
//	if (r && iter == end)
//	{
//		//			std::cout << "-------------------------\n";
//		//			std::cout << "Parsing succeeded\n";
//		//			//            std::cout << "result = " << result << std::endl;
//		//			std::cout << "-------------------------\n";
//	}
//	else
//	{
//		std::string rest(iter, end);
//		std::cout << "-------------------------\n";
//		std::cout << "Parsing failed\n";
//		std::cout << "stopped at: \": " << rest << "\"\n";
//		std::cout << "-------------------------\n";
//	}
//	return base_node;
	
//	bool r = phrase_parse(iter, end, func_parser,boost::spirit::ascii::space, parse_ret);

}


