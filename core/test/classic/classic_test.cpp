//This file is part of Bertini 2.
//
//classic_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//classic_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with classic_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
// classic_test.cpp:  main source file for the classic compatibility testing executable for Bertini2






//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Classic Compatibility Testing"
#include <boost/test/unit_test.hpp>





//#include "bertini.hpp"
//#include "classic.hpp"
//
//
//using namespace bertini::classic::parsing;
//using SplitFile = bertini::classic::SplitInputFile;



BOOST_AUTO_TEST_SUITE(classic_parsing_test)

//BOOST_AUTO_TEST_CASE(parser_works)
//{
//    std::string str = "CONFIG \n config info \n END \n INPUT \n input info \n END";
//    std::cout << str << std::endl;
//    
//    SplitFile split;
//    std::string::const_iterator iter = str.begin();
//    std::string::const_iterator end = str.end();
//    SplitFileInputConfig<std::string::const_iterator> S;
//    bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, split);
//    
//    BOOST_CHECK(s);
//    BOOST_CHECK(s && iter==end);
//}


BOOST_AUTO_TEST_SUITE_END()
