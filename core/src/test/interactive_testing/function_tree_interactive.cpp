//  function_tree_interactive.cpp
//  an interactive driver for the function and parser parts of bertini2.


#include <cmath>
#include <iostream>
#include <vector>


#include "function_tree.hpp"

#include <boost/spirit/include/qi.hpp>



namespace qi = boost::spirit::qi;



int InteractiveTest();




// set the values of the variables, by asking the user for their values.
void set_variable_values(std::vector<std::shared_ptr<Variable> > & variables);

// print variable values to the screen
void print_variable_values(const std::vector<std::shared_ptr<Variable> > & variables);

// get the names of the variables from the user.
std::vector<std::shared_ptr<Variable> > get_var_group_from_user(qi::symbols<char,int> *variable_names);


std::shared_ptr<Node> get_function_from_user(const std::vector<std::shared_ptr<Variable> > & variables, qi::symbols<char,int>  variable_names);






int main()
{
	return InteractiveTest();
	
}









int InteractiveTest()
{
	
	
	std::cout.precision(10);
	
	// 1. Read in the variables //
	std::cout << "/////////////////////////////////////////////////////////\n\n";
	std::cout << "\t\tPolynomial Playground\n\n";
	std::cout << "/////////////////////////////////////////////////////////\n\n";
	

	
	
	
	// 2. Main loop //

	std::vector<std::shared_ptr<Variable> > variables;

	
	qi::symbols<char,int> variable_names;
	
	std::shared_ptr<Node> F;
	
	
	bool cont = true;
	while (cont)
	{
		if (variables.size()){
			std::cout << "current variables: ";
			for (auto iter = variables.begin(); iter!=variables.end(); iter++) {
				std::cout << **iter << (iter==variables.end()-1 ? "\n" : ", ");
			}
		}
		else
			std::cout << "no variables declared\n";
		
		if (F){
			std::cout << "current function: " << *F << "\n";
		}
		else
			std::cout << "no function defined\n";
		
		std::cout << "\n   your options:  (case irrelevant)\n\n";
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
				variables = get_var_group_from_user(&variable_names);
				F.reset();
				break;
			}
			case 's':
			{
				set_variable_values(variables);
				if (F) {
					F->Reset();
				}
				break;
			}
			case 'e':
			{
				if (F){
					std::cout << "f = " << *F << std::endl;
					std::cout << "f(";
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
					F = get_function_from_user(variables,variable_names);
				}
				else
					std::cout << "no variables declared.  declare variables first." << std::endl;
				break;
			}
				
			case 'c':
			{
				F.reset();
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
	
	auto success = qi::phrase_parse(iter,end,var_parser.start(),boost::spirit::ascii::space);
	
	
	if (iter==end && success) {
		*variable_names = var_parser.variable();
		
		std::cout << "got 'em, thanks!" << std::endl;
	}
	else{
		std::cout << "parsing of the variables failed.\n";
		std::cout << "the unparsed part of the string is:\n\t" << std::string(iter,end) << std::endl;
	}
	
	
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
	
}

