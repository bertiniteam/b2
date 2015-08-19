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
//  ummer 2015

#include "bertini.hpp"
#include <string>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(classic_parsing)

BOOST_AUTO_TEST_CASE(split_apart_config_and_input_just_config_throws)
{
	std::string test_string = "CONFIG\n\ntracktype:1;\n\nEND;\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;

	bertini::classic::SplitInputFile config_and_input;

	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();

	BOOST_CHECK_THROW(phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input),std::runtime_error);

}

BOOST_AUTO_TEST_CASE(split_apart_config_and_input_just_input)
{
	std::string test_string = "INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;

	bertini::classic::SplitInputFile config_and_input;

	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);

	std::cout << config_and_input.Input() << "\n";

	std::size_t found_config = config_and_input.Input().find("variable_group x, y;");
	BOOST_CHECK(found_config!=std::string::npos);

	found_config = config_and_input.Input().find("INPUT");
	BOOST_CHECK(found_config==std::string::npos);

	found_config = config_and_input.Input().find("END;");
	BOOST_CHECK(found_config==std::string::npos);
}

BOOST_AUTO_TEST_CASE(split_apart_config_and_input_just_input_no_INPUT)
{
	std::string test_string = "\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;

	bertini::classic::SplitInputFile config_and_input;

	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);

	std::cout << config_and_input.Input() << "\n";

	std::size_t found_config = config_and_input.Input().find("variable_group x, y;");
	BOOST_CHECK(found_config!=std::string::npos);

	found_config = config_and_input.Input().find("INPUT");
	BOOST_CHECK(found_config==std::string::npos);

	found_config = config_and_input.Input().find("END;");
	BOOST_CHECK(found_config==std::string::npos);
}

BOOST_AUTO_TEST_CASE(split_apart_config_and_input_have_both_check_config)
{
	std::string test_string = "CONFIG\n\ntracktype:1;\n\nEND;\n\nINPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;

	bertini::classic::SplitInputFile config_and_input;

	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);

	std::cout << "config:\n\n";
	std::cout << config_and_input.Config() << "\n";

	std::size_t found_config = config_and_input.Config().find("tracktype:1;");
	BOOST_CHECK(found_config!=std::string::npos);
}

BOOST_AUTO_TEST_CASE(split_apart_config_and_input_have_both_check_input)
{
	std::string test_string = "CONFIG\n\ntracktype:1;\n\nEND;\n\nINPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;\n\n";

	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;

	bertini::classic::SplitInputFile config_and_input;

	std::string::const_iterator iter = test_string.begin();
	std::string::const_iterator end = test_string.end();
	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);

	std::cout << "input:\n\n";
	std::cout << config_and_input.Input() << "\n";


	std::size_t found_input = config_and_input.Input().find("variable_group x, y;");
	BOOST_CHECK(found_input!=std::string::npos);

}

BOOST_AUTO_TEST_SUITE_END()



