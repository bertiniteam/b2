//This file is part of Bertini 2.
//
//classic_parsing_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//classic_parsing_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with classic_parsing_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

//  parsing.cpp
//
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#include "bertini.hpp"
#include <bertini2/io/parsing/classic_utilities.hpp>
#include <string>
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(classic_parsing)


BOOST_AUTO_TEST_CASE(config_end_input_end_15) // 15
{
	std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\nEND; between INPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());

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

//	bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//	bertini::classic::SplitInputFile config_and_input;
//	std::string::const_iterator iter = test_string.begin();
//	std::string::const_iterator end = test_string.end();
//	phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
	auto config = config_and_input.Config();
	auto input = config_and_input.Input();

    BOOST_CHECK(config_and_input.Readable());


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

BOOST_AUTO_TEST_CASE(config_end_9) // 9
{
    std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n between \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;  efgh";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    
    BOOST_CHECK(config.find("abcd")==std::string::npos);
    BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
    BOOST_CHECK(config.find("variable_group x, y;")!=std::string::npos);
    BOOST_CHECK(config.find("f = x^2 + y^2 - 1;")!=std::string::npos);
    BOOST_CHECK(config.find("between")!=std::string::npos);
    
    
    BOOST_CHECK(config.find("CONFIG")==std::string::npos);
    BOOST_CHECK(config.find("INPUT;")==std::string::npos);
    BOOST_CHECK(config.find("END;")==std::string::npos);
    
    
}


BOOST_AUTO_TEST_CASE(config__8) // 8
{
    std::string test_string = "abcd CONFIG\n\ntracktype: 1;\n\n between \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n  efgh";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    
    BOOST_CHECK(!config_and_input.Readable());
    
    
}


BOOST_AUTO_TEST_CASE(end_input_end_7) // 7
{
    std::string test_string = "abcd tracktype: 1;\n\n END; between INPUT \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n END; efgh";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    
    BOOST_CHECK(!config_and_input.Readable());
    
    
}

BOOST_AUTO_TEST_CASE(end_input_6) // 6
{
    std::string test_string = "abcd tracktype: 1;\n\n END; between INPUT \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\n efgh";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    
    BOOST_CHECK(!config_and_input.Readable());
    
    
}


BOOST_AUTO_TEST_CASE(end___end_5) // 5
{
    std::string test_string = "abcd tracktype: 1;\n\n END; between \n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\nEND; \n efgh";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> parser;
//    bertini::classic::SplitInputFile config_and_input;
//    std::string::const_iterator iter = test_string.begin();
//    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);
	
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    
    BOOST_CHECK(!config_and_input.Readable());
    
    
}


BOOST_AUTO_TEST_CASE(uncomment)
{
    std::string test_string = "%Title of file\n  \n tracktype: 1;  %comment about setting\n %  More full comments\n %Another line of comments\n trackit: 12;\n %commentsetting: 4; \n  %%%%%%%%%%%%%%%%%END of Settings%%%%%%%%%%%%%%\n";
    
    bertini::parsing::classic::CommentStripper<std::string::const_iterator> parser;
	bertini::parsing::classic::SplitInputFile config_and_input;
    std::string::const_iterator iter = test_string.begin();
    std::string::const_iterator end = test_string.end();
    
    std::string test_out = "";
    bool s = phrase_parse(iter, end, parser, boost::spirit::ascii::space, test_out);
    
    std::string rest(iter, end);
    
    
    BOOST_CHECK(s && iter==end);
    
    BOOST_CHECK(test_out.find("%")==std::string::npos);
    BOOST_CHECK(test_out.find("Title of file")==std::string::npos);
    BOOST_CHECK(test_out.find("comment about setting")==std::string::npos);
    BOOST_CHECK(test_out.find("Another line of comments")==std::string::npos);
    BOOST_CHECK(test_out.find("commentsetting: 4;")==std::string::npos);
    BOOST_CHECK(test_out.find("END of Settings")==std::string::npos);
    
    BOOST_CHECK(test_out.find("tracktype: 1;")!=std::string::npos);
    BOOST_CHECK(test_out.find("trackit: 12;")!=std::string::npos);
}


BOOST_AUTO_TEST_CASE(test_split_and_uncomment)
{
    std::string test_string = "%Title of file\n CONFIG \n tracktype: 1;  %comment about setting\n %  More full comments\n %Another line of comments\n trackit: 12;\n %commentsetting: 4; \n  %%%%%%%%%%%%%%%%%END of Settings%%%%%%%%%%%%%%\n END; \n stuff %more comments\n INPUT\n %Beginning comments\n variable_group x,y; %variables\n % Parameters \n parameter t; \n function f\n %Polynomials \n f = x^2 + y;\n %End of INput\n END; stuff end";
    
//    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> split_parser;
    bertini::parsing::classic::CommentStripper<std::string::const_iterator> comment_parser;
//    bertini::classic::SplitInputFile config_and_input;
    std::string::const_iterator iter = test_string.begin();
    std::string::const_iterator end = test_string.end();
//    phrase_parse(iter, end, split_parser, boost::spirit::ascii::space, config_and_input);
	
	auto config_and_input = bertini::parsing::classic::ParseInputFile(test_string);

    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    std::string test_out = "";
    iter = config.begin();
    end = config.end();
    phrase_parse(iter, end, comment_parser, boost::spirit::ascii::space, test_out);
    config_and_input.SetConfig(test_out);
    
    test_out = "";
    iter = input.begin();
    end = input.end();
    phrase_parse(iter, end, comment_parser, boost::spirit::ascii::space, test_out);
    config_and_input.SetInput(test_out);
    config = config_and_input.Config();
    input = config_and_input.Input();
    
    
    std::string rest(iter, end);
    
    
    
    BOOST_CHECK(config.find("%")==std::string::npos);
    BOOST_CHECK(config.find("Title of file")==std::string::npos);
    BOOST_CHECK(config.find("comment about setting")==std::string::npos);
    BOOST_CHECK(input.find("%")==std::string::npos);
    BOOST_CHECK(input.find("Title of file")==std::string::npos);
    BOOST_CHECK(input.find("comment about setting")==std::string::npos);
    
    BOOST_CHECK(config.find("tracktype: 1;")!=std::string::npos);
    BOOST_CHECK(config.find("trackit: 12;")!=std::string::npos);

    BOOST_CHECK(input.find("variable_group x,y;")!=std::string::npos);
    BOOST_CHECK(input.find("f = x^2 + y;")!=std::string::npos);
    
}


BOOST_AUTO_TEST_SUITE_END()



