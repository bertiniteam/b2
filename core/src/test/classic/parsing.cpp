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


BOOST_AUTO_TEST_CASE(config_end_input_end)
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND; between INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

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
	BOOST_CHECK(input.find("efgh")==std::string::npos);
	BOOST_CHECK(input.find("between")==std::string::npos);

	BOOST_CHECK(config.find("CONFIG")==std::string::npos);
	BOOST_CHECK(config.find("INPUT;")==std::string::npos);
	BOOST_CHECK(config.find("END;")==std::string::npos);
	
	BOOST_CHECK(input.find("INPUT")==std::string::npos);
	BOOST_CHECK(input.find("CONFIG")==std::string::npos);
	BOOST_CHECK(input.find("END;")==std::string::npos);

}





BOOST_AUTO_TEST_CASE(config_input_end)
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

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


BOOST_AUTO_TEST_CASE(config_input)
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


BOOST_AUTO_TEST_CASE(config_end_end)
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



BOOST_AUTO_TEST_CASE(config_end__no_input_markers)
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








BOOST_AUTO_TEST_CASE(input_end)
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


BOOST_AUTO_TEST_CASE(input)
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

BOOST_AUTO_TEST_CASE(end)
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

BOOST_AUTO_TEST_SUITE_END()



