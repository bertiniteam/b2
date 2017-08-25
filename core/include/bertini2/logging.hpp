//This file is part of Bertini 2.
//
//logging.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//logging.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with logging.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file logging.hpp 

\brief Logging in Bertini using Boost.Log
*/


#ifndef BERTINI_LOGGING_HPP
#define BERTINI_LOGGING_HPP 

 
#define BOOST_LOG_DYN_LINK 1

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace bertini
{


	namespace logging = boost::log;
	namespace src = boost::log::sources;
	namespace sinks = boost::log::sinks;
	namespace keywords = boost::log::keywords;
	
	// the following is adapted from https://stackoverflow.com/questions/11421432/
	// question answered by user James Adkison, asked by Adi, edited by James McNellis.
	// the adaptation is the replacement of std::ostream with the blos type.  why the unmodified code still fails
	// for blos types is a mystery.
	using blos = boost::log::record_ostream;
	template<typename T>
	blos& operator<<(typename std::enable_if<std::is_enum<T>::value, blos>::type& stream, const T& e)
	{
		return stream << static_cast<typename std::underlying_type<T>::type>(e);
	}


	struct LoggingInit
	{
		
		// trivial logger-provided severity levels are  
		//
		//  trace, debug, info, warning, error, fatal
		LoggingInit(logging::trivial::severity_level desired_level = logging::trivial::severity_level::trace, unsigned desired_rotation_size = 10*1024*1024)
		{
			logging::add_file_log
			(
			    keywords::file_name = "bertini_%N.log",
			    keywords::rotation_size = desired_rotation_size,
			    keywords::format = "%Message%" //[%TimeStamp%]: 
			);

			logging::core::get()->set_filter
			(
			    logging::trivial::severity >= desired_level
			);

			BOOST_LOG_TRIVIAL(trace) << "initialized logging";

		}

		~LoggingInit(){}
	    
	};

	using severity_level = logging::trivial::severity_level;

} // re: namespace bertini



#endif
