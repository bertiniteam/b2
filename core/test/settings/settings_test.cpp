//This file is part of Bertini 2.0.
//
//function_tree_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_tree_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_tree_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  settings_test.cpp
//
//  Created by Collins, James B. on 8/27/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Settings Testing"


#include <iostream>

#include <cstdlib>
#include <cmath>

#include "bertini2/bertini.hpp"
#include "bertini2/settings/configIni_parse.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

using mpfr = bertini::mpfr;

using namespace boost::program_options;
int tracktype_id;
void SetTrackType(int type)
{
    tracktype_id = type;
}



BOOST_AUTO_TEST_SUITE(config_settings)


BOOST_AUTO_TEST_CASE(parse_store_set_settings)
{
    std::stringstream configfile("TrackType = 3  \n MPType=2\n Precision  =40");
    
    options_description options{"test settings"};
    options_description track_type{"Track Type Settings"};
    options_description MP{"Multiple Precision Settings"};
    
    track_type.add_options()
    ("TrackType", value<int>()->default_value(0)->notifier(&SetTrackType), "Defines the type of tracking that Bertini will perform");
    
    mpfr z;
    MP.add_options()
    ("MPType", value<int>()->default_value(2), "Defines the type of multiple precision used by Bertini")
    ("Precision", value<int>()->default_value(96)->notifier(boost::bind(&mpfr::precision,&z, _1)), "The precision to use in calculations");
    
    options.add(track_type);
    options.add(MP);
    
    variables_map vm;
    store(parse_config_file(configfile, options), vm);
    notify(vm);

    BOOST_CHECK(tracktype_id==3);
    BOOST_CHECK(vm["MPType"].as<int>() == 2);
    BOOST_CHECK_EQUAL(z.precision(),40);

}


BOOST_AUTO_TEST_CASE(parsing_to_ini)
{
    std::string test_string = "tracktype :1 \n MPType:-1";
    
    bertini::settings::parsing::ConfigToIni<std::string::const_iterator> parser;
    std::string::const_iterator iter = test_string.begin();
    std::string::const_iterator end = test_string.end();
    
    std::string str_out;
    phrase_parse(iter, end, parser, boost::spirit::ascii::space, str_out);
    std::cout << "d " << str_out << std::endl;
    BOOST_CHECK(str_out.find("tracktype=1")!=std::string::npos);
    BOOST_CHECK(str_out.find("mptype=-1")!=std::string::npos);

}


BOOST_AUTO_TEST_CASE(full_config_to_vm)
{
    std::string test_string = "%Title of file\n CONFIG \n TrackType: 3;  %comment about setting\n %  More full comments\n %Another line of comments\n PRECISION:40 \n  %commentsetting: 4; \n  %%%%%%%%%%%%%%%%%END of Settings%%%%%%%%%%%%%%\n END; \n stuff %more comments\n INPUT\n %Beginning comments\n variable_group x,y; %variables\n % Parameters \n parameter t; \n function f\n %Polynomials \n f = x^2 + y;\n %End of INput\n END; stuff end";
    
    bertini::classic::parsing::SplitFileInputConfig<std::string::const_iterator> split_parser;
    bertini::classic::parsing::CommentStripper<std::string::const_iterator> comment_parser;
    bertini::settings::parsing::ConfigToIni<std::string::const_iterator> config_ini_parser;
    bertini::classic::SplitInputFile config_and_input;
    std::string::const_iterator iter = test_string.begin();
    std::string::const_iterator end = test_string.end();
    phrase_parse(iter, end, split_parser, boost::spirit::ascii::space, config_and_input);
    auto config = config_and_input.Config();
    auto input = config_and_input.Input();
    
    std::string test_out = "";
    iter = config.begin();
    end = config.end();
    phrase_parse(iter, end, comment_parser, boost::spirit::ascii::space, test_out);
    config_and_input.SetConfig(test_out);
    config = config_and_input.Config();
    
    test_out = "";
    iter = config.begin();
    end = config.end();
    phrase_parse(iter, end, config_ini_parser, boost::spirit::ascii::space, test_out);
    config_and_input.SetConfig(test_out);
//    std::string rest(iter,end);
//    std::cout << "rest = " << test_out << std::endl;
    
    std::stringstream configfile(config_and_input.Config());
    
    options_description options{"test settings"};
    options_description track_type{"Track Type Settings"};
    options_description MP{"Multiple Precision Settings"};
    
    track_type.add_options()
    ("tracktype", value<int>()->default_value(0)->notifier(&SetTrackType), "Defines the type of tracking that Bertini will perform");
    
    mpfr z;
    MP.add_options()
    ("mptype", value<int>()->default_value(2), "Defines the type of multiple precision used by Bertini")
    ("precision", value<int>()->default_value(96)->notifier(boost::bind(&mpfr::precision,&z, _1)), "The precision to use in calculations");
    
    options.add(track_type);
    options.add(MP);
    
    variables_map vm;
    store(parse_config_file(configfile, options), vm);
    notify(vm);
    
    BOOST_CHECK(tracktype_id==3);
    BOOST_CHECK(vm["mptype"].as<int>() == 2);
    BOOST_CHECK_EQUAL(z.precision(),40);
    
}



BOOST_AUTO_TEST_SUITE_END()
