//This file is part of Bertini 2.
//
//bertini2/blackbox/algorithm_builder.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/algorithm_builder.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/algorithm_builder.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// individual authors of this file include 
//
// Dani Brake, university of notre dame

/**
\file bertini2/blackbox/algorithm_builder.hpp 

\brief A type which determines which algorithms etc, to run, and sets them up.
*/



#pragma once

#include "bertini2/io/file_utilities.hpp"

namespace bertini {
namespace blackbox {


/**
\brief A de-serializer, whose responsibility it is to construct runnable objects based on input files.

This class is inspired by the SimBuilder class from Hythem Sidky's SAPHRON package for molecular dynamics.
*/
class AlgoBuilder
{
	
public:
	AlgoBuilder() = default;

	/**
	\brief Method for constructing an algorithm from a bertini classic input file
	*/
	int ClassicBuild(boost::filesystem::path const& input_file);

	/**
	\brief Returns a non-owning pointer to the built algorithm
	*/
	AnyAlgorithm* GetAlg()
	{
		return alg_.get();
	}

private:
	std::unique_ptr<AnyAlgorithm> alg_;
};



} // blackbox
} //namespace bertini

