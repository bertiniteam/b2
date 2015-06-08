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
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015



#include <boost/test/unit_test.hpp>



#include "system.hpp"
#include "system_parsing.hpp"

using System = bertini::System;


BOOST_AUTO_TEST_SUITE(system_class)


BOOST_AUTO_TEST_CASE(system_make_a_system_at_all)
{
	System S;
}

BOOST_AUTO_TEST_CASE(system_create_parser)
{
	System sys;
	std::string str = "variable_group x, y, z;\nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";
	
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	
	
	bertini::SystemParser<std::string::const_iterator> S;
	
	
	bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);
	
	BOOST_CHECK(s && iter==end);
	
}




BOOST_AUTO_TEST_SUITE_END()




