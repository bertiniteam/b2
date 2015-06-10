//This file is part of Bertini 2.0.
//
//bertini.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
// bertini.cpp:  main source file for computational core of Bertini2



#include "bertini.hpp"


int Node::tabcount = 0;



int main(int argc, char** argv)
{
	std::cout << "Bertini, version 2.\n\nDeveloped by B-Team.\n\nThis program is currently under development,\nand doesn't do anything useful yet." << std::endl;

	
	bertini::System sys;
	std::string str = "variable_group x, y, z; function f1, f2; pathvariable t ; parameter p, q;";
	
	
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	
	
	bertini::SystemParser<std::string::const_iterator> S;
	
	
	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
	
	if (!s || iter!=end){
		std::cout << "parsing failed.\n";
		std::cout << "the unparsed string:\n" << std::string(iter,end) << "\n";
	}
	else{
		std::cout << "parsing of the entire string was successful\n";
	}
	
	
	std::cout << sys << std::endl;
	return 0;
}








