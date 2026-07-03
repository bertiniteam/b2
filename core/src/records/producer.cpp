//This file is part of Bertini 2.
//
//producer.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//producer.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with producer.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file producer.cpp

\brief Producer identification for run records.  The version comes from the
BERTINI2_VERSION_FULL compile definition (the VERSION file); the commit from a
build-time generated header (cmake/GenerateProducerInfo.cmake), so it tracks the
actual source state, not the last cmake configure.
*/

#include "bertini2/records/producer.hpp"

#include "bertini2/records/producer_info.generated.hpp"

namespace bertini {
namespace records {

std::string ProducerVersion()
{
#ifdef BERTINI2_VERSION_FULL
	return BERTINI2_VERSION_FULL;
#else
	return "unknown";
#endif
}

std::string ProducerCommit()
{
	return BERTINI2_GIT_COMMIT;
}

boost::json::object ProducerInfo()
{
	return {{"name", "bertini2"}, {"version", ProducerVersion()}, {"commit", ProducerCommit()}};
}

} // namespace records
} // namespace bertini
