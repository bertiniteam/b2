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

#include <iostream>

#include <cstdlib>
#include <cmath>

#include "bertini.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/program_options.hpp>

using namespace boost::program_options;
int tracktype_id;
void SetTrackType(int type)
{
    tracktype_id = type;
}


BOOST_AUTO_TEST_SUITE(config_settings)


BOOST_AUTO_TEST_CASE(initial_test)
{
    
    
    
    std::stringstream configfile("   \n");
    
    options_description options{"test settings"};
    options_description track_type{"Track Type Settings"};
    options_description MP{"Multiple Precision Settings"};
    
    track_type.add_options()
    ("Track Type", value<int>()->default_value(0)->notifier(&SetTrackType), "Defines the type of tracking that Bertini will perform");
    
    mpfr z;
    MP.add_options()
    ("MPType", value<int>()->default_value(2), "Defines the type of multiple precision used by Bertini")
    ("Precision", value<int>()->default_value(96)->notifier(&z.precision), "The precision to use in calculations");
    
    variables_map vm;
    store(parse_config_file(configfile, track_type), vm);
    notify(vm);

    BOOST_CHECK_EQUAL(tracktype_id,0);

}





BOOST_AUTO_TEST_SUITE_END()
