//This file is part of Bertini 2.
//
//bertini2/io/splash.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/splash.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/splash.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/io/splash.hpp 

\brief Provides the splash screens for bertini2.
*/

#pragma once

#include "bertini2/config.h"

#include "boost/version.hpp"

#include <sstream>

namespace bertini{

inline
std::string LicenseInfo()
{
	std::stringstream ss;
	ss << "Bertini 2 is available under the GPL 3 license,\nwith additional terms as permitted under Section 7.\n\nPlease see the full text of the licenses for Bertini and its dependencies in the source code at b2/licenses.\n\n";
	ss << "Bertini 2 depends on the following software and libraries:"
	"\n"
	"==\n"
	"Eigen\n"
	"--\n"
	"Eigen is available from http://eigen.tuxfamily.org/index.php?title=Main_Page\n"
	"\n"
	"Eigen's license as of April 12, 2016 is MPL2, with some LGPL scattered in.\nFor more information, please consult\n"
	"http://eigen.tuxfamily.org/index.php?title=FAQ#Licensing\n"
	"\n"
	"The Mozilla Public License 2.0 has been included at b2/licenses/\n"
	"\n"
	"==\n"
	"GMP\n"
	"--\n"
	"The GNU Multiple Precision Library is available from https://gmplib.org/.\n"
	"\n"
	"Since version 6, GMP is distributed under the dual licenses, GNU LGPL v3 and GNU GPL v2.\n"
	"\n"
	"The GNU LGPL license has been included in Bertini 2 at b2/licenses/\n"
	"\n"
	"==\n"
	"MPFR\n"
	"--\n"
	"The GNU MPFR library is available from http://www.mpfr.org/.\n"
	"\n"
	"MPFR is free.\nIt is distributed under the GNU Lesser General Public License (GNU Lesser GPL),\nversion 3 or later (2.1 or later for MPFR versions until 2.4.x).\n"
	"\n"
	"==\n"
	"Boost\n"
	"--\n"
	"The Boost C++ libraries are available from http://www.boost.org/.\n"
	"\n"
	"Boost is available under the Boost License 1.0.\nPlease see http://www.boost.org/users/license.html.\nA copy of the BPL has been included in b2/licenses.\n";

	return ss.str();
}

inline 
std::string SourceURL()
{
	return PACKAGE_URL;
}

inline
std::string WikiURL()
{
	return "https://github.com/bertiniteam/b2/wiki";
}

inline
std::string Version()
{
	return PACKAGE_VERSION;
}

inline 
std::string Owners()
{
	std::stringstream ss;
	ss << "D.J. Bates, D.A. Brake, J.D. Hauenstein,\nA.J. Sommese, C.W. Wampler";
	return ss.str();
}

inline 
std::string Authors()
{
	std::stringstream ss;
	ss << "D. Brake, J. Collins, T. Hodges";
	return ss.str();
}

inline
std::string SplashScreen()
{
	std::stringstream ss;

	ss << "    Bertini(TM) 2\n\n";
	ss << "  The Bertini Trademark is owned by\n" << Owners() << "\n\n";
	ss << "  The code is primarily authored by\n" << Authors() << "\n\n";
	ss << "  Source available online at\n" << SourceURL() << "\n\n";
	ss << "  Wiki online at\n" << WikiURL() << "\n\n";
	ss << "  This is version\n" << Version() << "\n";
	ss << "  Bertini2 is GPL3 Free/Libre Open Source Software, please contribute!\n\n";
	return ss.str();
}


inline 
std::string GenericHelp()
{
	return "";
}


inline 
std::string BoostHeaderVersion()
{
	std::stringstream ss;
	ss << BOOST_VERSION; 
	return ss.str();
}



inline 
std::string DependencyVersions()
{
	std::stringstream ss;
	ss << "Compiled against Boost headers " << BoostHeaderVersion() << "\n\n";
	return ss.str();
}


}// re: namespace bertini


