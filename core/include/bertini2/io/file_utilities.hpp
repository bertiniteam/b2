//This file is part of Bertini 2.
//
//bertini2/io/file_utilities.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/file_utilities.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/file_utilities.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/io/file_utilities.hpp 

\brief Provides the file_utilities screens for bertini2.
*/

#pragma once
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

namespace bertini{

	namespace fs = boost::filesystem;

	using Path = fs::path;

	using ifstream = std::ifstream;
	/**
	\brief Try to open a file, and throw if it doesn't exist, or is a directory.
	*/
	inline
	void OpenInFileThrowIfFail(ifstream & in, 
	                           Path const& input_path)
	{
		using namespace fs;
		try{
			if (exists(input_path))
			{
				if (is_directory(input_path))
				{
					std::stringstream err_msg;
					err_msg << "attempting to open file but is a directory, of name '" << input_path.string() << "'";
					throw std::runtime_error(err_msg.str());
				}
				else
				{
					in.open(input_path.c_str());
					if (!in.is_open())
					{
						std::stringstream err_msg;
						err_msg << "file '" << input_path.string() << "' hypothetically exists, but failed to open correctly";
						throw std::runtime_error(err_msg.str());
					}
				}
			} 
			else
			{
				std::stringstream err_msg;
				err_msg << "attempting to open file which doesn't exist, of name '" << input_path.string() << "'";
				throw std::runtime_error(err_msg.str());
			}
		}
		catch (const filesystem_error& ex)
		{
			std::stringstream err_msg;
			err_msg << "boost::filesystem throw while attempting to open file '" << input_path.string() << "': '" << ex.what() << "'";
			throw std::runtime_error(err_msg.str());
		}
	}

	/**
	\brief Read an entire file into a string.

	\return The string, now contaning the file.
	\param input_path The path to the file.
	*/
	inline
	std::string FileToString(Path const& input_path)
	{
		ifstream infile;
		OpenInFileThrowIfFail(infile, input_path);
		return std::string ( std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>() );
	}
}

