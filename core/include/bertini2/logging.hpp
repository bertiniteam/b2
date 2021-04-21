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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

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

namespace bertini{
namespace logging{

	

	namespace blog = boost::log;

	using severity_level = blog::trivial::severity_level;
	namespace src = blog::sources;
	namespace sinks = blog::sinks;
	namespace keywords = blog::keywords;
	

	// the following is adapted from https://stackoverflow.com/questions/11421432/
	// question answered by user James Adkison, asked by Adi, edited by James McNellis.
	// the adaptation is the replacement of std::ostream with the blros type.  why the unmodified code still fails
	// for blros types is a mystery.
	using blros = blog::record_ostream;
	template<typename T>
	blros& operator<<(typename std::enable_if<std::is_enum<T>::value, blros>::type& stream, const T& e)
	{
		return stream << static_cast<typename std::underlying_type<T>::type>(e);
	}


	/**
	\class Logging

	Provided as an interface to the underlying logging library.

	Highlight functions:
	
	* Init -- a "call-it-once" kinda function
	* SetFilter
	* AddFile

	I have no idea how to remove a file, once you have done AddFile.  If this is something you need, contact silviana amethyst, and ask her to provide such a function.  She practices YAGNI, and she hadn't NI yet.

	There is Init, with all defaults, so you should totally call it to initialize all logging facilities for Bertini2.  Failure to do so produces pure screen output.

	//[%TimeStamp%]: 
	*/
	struct Logging
	{

		static 
		void Init(std::string const& name_pattern = "bertini_%N.log", 
				std::string const& format = "%Message%", 
				unsigned rotation_size = 10*1024*1024,
				severity_level const& new_level = severity_level::error)
		{
			AddFile(name_pattern, format, rotation_size);
			SetLevel(new_level);

			BOOST_LOG_TRIVIAL(info) << "initialized logging";
		}


		static
		void AddFile(std::string const& name_pattern, std::string const& format, unsigned rotation_size, bool auto_flush = true)
		{
			blog::add_file_log
			(
			    keywords::file_name = name_pattern,
			    keywords::rotation_size = rotation_size,
			    keywords::format = format,
			    keywords::auto_flush = auto_flush 
			);
		}


		/**
		 trivial logger-provided severity levels are  

		 trace, debug, info, warning, error, fatal
		*/
		static
		void SetLevel(severity_level const& new_level)
		{
			blog::core::get()->set_filter
			(
			    blog::trivial::severity >= new_level
			);
		}
		
	    
	};

	

} // namespace logging
} // re: namespace bertini



#endif
