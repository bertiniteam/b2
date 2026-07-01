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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

//  parsing.cpp


#include "bertini2/bertini.hpp"
#include <bertini2/io/parsing/classic_utilities.hpp>
#include <bertini2/io/parsing/system_parsers.hpp>
#include <bertini2/io/classic_writer.hpp>
#include <string>
#include <set>
#include <fstream>
#include <boost/filesystem.hpp>
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


// The classic WRITER is the inverse of the parser: emitting a system to classic syntax and parsing
// it back must reconstruct an equivalent system (so a Bertini 1 run sees the same problem).
BOOST_AUTO_TEST_CASE(classic_writer_round_trips_a_system)
{
	using namespace bertini;
	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");
	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x*x + y*y - 1);
	sys.AddFunction(x - y);

	System reparsed{ classic::SystemToClassic(sys) };
	BOOST_CHECK_EQUAL(reparsed.NumNaturalFunctions(), sys.NumNaturalFunctions());

	// identical values at a generic point -> the emitted classic text is a faithful round-trip
	Vec<complex_dbl> pt(2); pt << complex_dbl(0.3, 0.7), complex_dbl(-0.4, 0.2);
	auto a = sys.Eval(pt);
	auto b = reparsed.Eval(pt);
	BOOST_REQUIRE_EQUAL(a.size(), b.size());
	for (Eigen::Index i = 0; i < a.size(); ++i)
		BOOST_CHECK_SMALL(abs(a(i) - b(i)), 1e-12);
}


// Regression: a classic input FILE with '%' comments must parse.  The blackbox reads files via
// the Path overload of SplitIntoConfigAndInput, which used to split the raw text WITHOUT running
// the CommentStripper first -- so a '%' comment (Bertini 1's comment marker) survived into the
// input section and the system parser choked ("did not consume entire input").  Comments appear
// here both on their own line and trailing real declarations, in both CONFIG and INPUT.
BOOST_AUTO_TEST_CASE(file_with_percent_comments_parses)
{
	namespace fs = boost::filesystem;
	auto path = fs::temp_directory_path() / fs::unique_path("b2_comment_%%%%-%%%%.b2");

	{
		std::ofstream out(path.string());
		out <<
			"CONFIG\n"
			"% a full-line comment in config\n"
			"tracktype: 0;   % trailing comment after a real setting\n"
			"END;\n"
			"INPUT\n"
			"% two groups => multihomogeneous; this comment must be stripped\n"
			"variable_group x;\n"
			"variable_group y;   % trailing comment on a declaration\n"
			"function f1, f2;\n"
			"f1 = x*y - 1;\n"
			"f2 = x + y - 3;   % and another\n"
			"END;\n";
	}

	std::string config, input;
	bertini::parsing::classic::SplitIntoConfigAndInput(config, input, path);
	fs::remove(path);

	// the '%' comment text is gone from both sections...
	BOOST_CHECK(config.find('%') == std::string::npos);
	BOOST_CHECK(input.find('%')  == std::string::npos);
	// ...while the real declarations survive...
	BOOST_CHECK(input.find("variable_group x;") != std::string::npos);
	BOOST_CHECK(input.find("f1 = x*y - 1;")     != std::string::npos);

	// ...and the system parser actually consumes the comment-free input.
	bertini::System sys;
	auto iter = input.begin();
	auto end  = input.end();
	bertini::parsing::classic::parse(iter, end, sys);
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2u);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
}


// ---- Unicode (UTF-8) identifier support ----
// The classic grammar accepts Unicode *letters* as identifiers (Ω, α, CJK, ...),
// storing names as their raw UTF-8 bytes.  UTF-8 is spelled with byte escapes so
// the source stays plain ASCII (portable across compilers), and inputs are built
// by std::string concatenation so a hex escape never swallows a following digit.
// See io/parsing/unicode_ident.hpp.
namespace {
	const std::string kOmega = "\xCE\xA9";         // U+03A9 GREEK CAPITAL LETTER OMEGA
	const std::string kAlpha = "\xCE\xB1";         // U+03B1 GREEK SMALL LETTER ALPHA
	const std::string kCJK   = "\xE4\xB8\xAD";     // U+4E2D
	const std::string kParty = "\xF0\x9F\x8E\x89"; // U+1F389 PARTY POPPER (single-code-point emoji)
	// 👍🏽 = THUMBS UP (U+1F44D) + skin-tone modifier (U+1F3FD): a two-code-point emoji.
	const std::string kThumbsToned = "\xF0\x9F\x91\x8D" "\xF0\x9F\x8F\xBD";
	// 👩‍👩‍👧 = WOMAN + ZWJ + WOMAN + ZWJ + GIRL: a five-code-point ZWJ sequence.
	const std::string kFamily = "\xF0\x9F\x91\xA9" "\xE2\x80\x8D" "\xF0\x9F\x91\xA9"
	                            "\xE2\x80\x8D" "\xF0\x9F\x91\xA7";
}

BOOST_AUTO_TEST_CASE(unicode_variable_group_omega_parses)
{
	std::string input = "variable_group " + kOmega + ", " + kAlpha + ";\n"
	                    "function f;\n"
	                    "f = " + kOmega + "^2 + " + kAlpha + "^2 - 1;\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 1u);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	auto names = sys.VariableNameSet();
	BOOST_CHECK(names.count(kOmega) == 1);
	BOOST_CHECK(names.count(kAlpha) == 1);
}

BOOST_AUTO_TEST_CASE(unicode_name_roundtrips_utf8)
{
	using namespace bertini;
	std::string input = "variable_group " + kOmega + ";\n"
	                    "function f;\n"
	                    "f = " + kOmega + "^2 - 2;\n";
	System sys{ input };
	std::string emitted = classic::SystemToClassic(sys);   // the classic writer emits UTF-8
	BOOST_CHECK(emitted.find(kOmega) != std::string::npos);
}

BOOST_AUTO_TEST_CASE(unicode_cjk_variable_parses)
{
	std::string input = "variable_group " + kCJK + ", y;\n"
	                    "function f;\n"
	                    "f = " + kCJK + " + y;\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	BOOST_CHECK(sys.VariableNameSet().count(kCJK) == 1);
}

BOOST_AUTO_TEST_CASE(unicode_mixed_ascii_unicode_identifier)
{
	// A single identifier mixing ASCII and Unicode letters/digits.
	std::string ident = "x" + kOmega + "1";
	std::string input = "variable_group " + ident + ", y;\n"
	                    "function f;\n"
	                    "f = " + ident + " + y;\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	BOOST_CHECK(sys.VariableNameSet().count(ident) == 1);
}

BOOST_AUTO_TEST_CASE(unicode_symbol_prefix_not_greedy)
{
	// Ω is declared but Ωα is not; the boundary guard stops the known symbol Ω from
	// matching a prefix of Ωα, so referencing Ωα must fail to parse rather than
	// silently becoming Ω (leaving α dangling).
	std::string input = "variable_group " + kOmega + ";\n"
	                    "function f;\n"
	                    "f = " + kOmega + kAlpha + ";\n";
	BOOST_CHECK_THROW(bertini::System sys{ input }, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(utf8_bom_is_stripped)
{
	// A leading UTF-8 BOM must be dropped by System(std::string) so it is not seen
	// as a stray leading character.
	std::string input = std::string("\xEF\xBB\xBF") +
	                    "variable_group x, y;\n"
	                    "function f;\n"
	                    "f = x^2 + y^2 - 1;\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 1u);
}

BOOST_AUTO_TEST_CASE(emoji_variable_parses)
{
	// Single-code-point emoji are valid identifiers.
	std::string input = "variable_group " + kParty + ", y;\n"
	                    "function f;\n"
	                    "f = " + kParty + " + y;\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	BOOST_CHECK(sys.VariableNameSet().count(kParty) == 1);
}

BOOST_AUTO_TEST_CASE(multipoint_emoji_variable_parses)
{
	// A multi-code-point emoji (skin-toned, and a ZWJ sequence) reads as ONE
	// contiguous identifier -- the whole byte sequence is the variable's name.
	std::string input = "variable_group " + kThumbsToned + ", " + kFamily + ";\n"
	                    "function f;\n"
	                    "f = " + kThumbsToned + " + " + kFamily + ";\n";
	bertini::System sys{ input };
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	auto names = sys.VariableNameSet();
	BOOST_CHECK(names.count(kThumbsToned) == 1);   // not split at the skin-tone modifier
	BOOST_CHECK(names.count(kFamily) == 1);        // not split at the ZWJs
}


BOOST_AUTO_TEST_SUITE_END()



