//This file is part of Bertini 2.0.
//
//parsing.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//parsing.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with parsing.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  parsing.cpp
//
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#include "bertini.hpp"
#include <string>
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(classic_parsing)


BOOST_AUTO_TEST_CASE(config_end_input_end_15) // 15
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND; between INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);
	BOOST_CHECK(config.find("between")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
    BOOST_CHECK(input.find("f = x^2 + y^2 - 1;")!=std::string::npos);
	BOOST_CHECK(input.find("efgh")==std::string::npos);
	BOOST_CHECK(input.find("between")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}


BOOST_AUTO_TEST_CASE(config_end_input_14) // 14
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND; between INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);
	BOOST_CHECK(config.find("between")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
	BOOST_CHECK(input.find("between")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}


BOOST_AUTO_TEST_CASE(config_end_end_13) // 13
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND;\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
	BOOST_CHECK(input.find("efgh")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}



BOOST_AUTO_TEST_CASE(config_end__no_input_markers_12) // 12
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND;\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}


BOOST_AUTO_TEST_CASE(config_input_end_11) // 11
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nINPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
	BOOST_CHECK(input.find("efgh")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}









BOOST_AUTO_TEST_CASE(config_input_10) // 10
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("tracktype: 1;")==std::string::npos);
	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}






//////////
//
//  TEST CASES WHERE CONFIG SHOULD BE EMPTY.
//
///////









BOOST_AUTO_TEST_CASE(input_end_3) // 3
{
	std::string test_string = "abcd INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
	BOOST_CHECK(input.find("efgh")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}


BOOST_AUTO_TEST_CASE(input_2) // 2
{
	std::string test_string = "abcd INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();


	BOOST_CHECK(config.find("abcd")==std::string::npos);
	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}

BOOST_AUTO_TEST_CASE(end_1) // 1
{
	std::string test_string = "\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();



	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);
	BOOST_CHECK(input.find("efgh")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}


BOOST_AUTO_TEST_CASE(no_markers_0) // 0
{
	std::string test_string = "\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
	bertini::classic::SplitInputFile config_and_input;
	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();



	BOOST_CHECK(config.find("variable_group x, y;")==std::string::npos);

	BOOST_CHECK(input.find("variable_group x, y;")!=std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}








////////////////
//
//  TEST CASES WHERE CONFIG IS DISCARDED OR THE INPUT FILE IS MALFORMED
//
//////////////////////

//BOOST_AUTO_TEST_CASE(config_end_9) // 9
//{
//    std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n between \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";
//    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
//    auto config = config_and_input.Config();
//    auto input = config_and_input.Input();
//    
//    
//    BOOST_CHECK(config.find("abcd")==std::string::npos);
//    BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("variable_group x, y;")!=std::string::npos);
//    BOOST_CHECK(config.find("f = x^2 + y^2 - 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("between")!=std::string::npos);
//    
//    
//    BOOST_CHECK(config.find("CONFIG")==std::string::npos);
//    BOOST_CHECK(config.find("INPUT;")==std::string::npos);
//    BOOST_CHECK(config.find("END;")==std::string::npos);
//    
//    
//    std::cout << input << std::endl;
//}
//
//
//BOOST_AUTO_TEST_CASE(config__8) // 8
//{
//    std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n between \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n  efgh";
//    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
//    auto config = config_and_input.Config();
//    auto input = config_and_input.Input();
//    
//    
//    BOOST_CHECK(config.find("abcd")==std::string::npos);
//    BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("variable_group x, y;")!=std::string::npos);
//    BOOST_CHECK(config.find("f = x^2 + y^2 - 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("efgh")!=std::string::npos);
//    
//    
//    BOOST_CHECK(config.find("CONFIG")==std::string::npos);
//    BOOST_CHECK(config.find("INPUT;")==std::string::npos);
//    BOOST_CHECK(config.find("END;")==std::string::npos);
//    
//    
//    std::cout << input << std::endl;
//}
//
//
//BOOST_AUTO_TEST_CASE(end_input_end_7) // 7
//{
//    std::string test_string = "abcd tracktype: 1;\n\n END; between INPUT \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n  efgh";
//    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
//    auto config = config_and_input.Config();
//    auto input = config_and_input.Input();
//    
//    
//    BOOST_CHECK(config.find("abcd")==std::string::npos);
//    BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("variable_group x, y;")!=std::string::npos);
//    BOOST_CHECK(config.find("f = x^2 + y^2 - 1;")!=std::string::npos);
//    BOOST_CHECK(config.find("efgh")!=std::string::npos);
//    
//    
//    BOOST_CHECK(config.find("CONFIG")==std::string::npos);
//    BOOST_CHECK(config.find("INPUT;")==std::string::npos);
//    BOOST_CHECK(config.find("END;")==std::string::npos);
//    
//    
//    std::cout << input << std::endl;
//}


BOOST_AUTO_TEST_SUITE_END()



