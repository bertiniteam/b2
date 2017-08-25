//This file is part of Bertini 2.
//
//bertini2/blackbox/algorithm_builder.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/algorithm_builder.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/algorithm_builder.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/blackbox/algorithm_builder.cpp 

\brief Provides the methods for building algorithms from files or streamable sources.
*/

namespace bertini{

namespace blackbox{

int AlgoBuilder::ClassicBuild(boost::filesystem::path const& input_file)
{
	std::string config, input;
	std::tie(config, input) = SplitIntoConfigAndInput(input_file);


	template<typename T>
	using AllConfs = blackbox::config::Configs::All<T>::type;
	auto results_double = bertini::parsing::classic::GetConfigSettings<double, AllConfs<double>>(config);
	auto results_mp = bertini::parsing::classic::GetConfigSettings<mpfr_float, AllConfs<mpfr_float>>(config);
	// now we have all the settings needed for any algorithm, huzzah



	return 0;
}



} // blackbox
} // bertini

