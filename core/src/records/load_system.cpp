//This file is part of Bertini 2.
//
//load_system.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//load_system.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with load_system.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file load_system.cpp

\brief Reloading an archived system: definition document -> canonical encoding ->
System, with the content digest as the acceptance test.
*/

#include "bertini2/records/load_system.hpp"

#include <boost/json.hpp>

#include <stdexcept>

namespace bertini {
namespace records {

std::string EncodingOfSystemDefinition(std::string const& definition_json)
{
    boost::json::value doc;
    try
    {
        doc = boost::json::parse(definition_json);
    }
    catch (std::exception const& e)
    {
        throw std::runtime_error(std::string("system definition: not JSON (") + e.what() + ")");
    }
    if (!doc.is_object())
        throw std::runtime_error("system definition: not a JSON object");
    auto const* encoding = doc.get_object().if_contains("encoding");
    if (!encoding || !encoding->is_string())
        throw std::runtime_error("system definition: no \"encoding\" string in the document");
    return std::string(encoding->get_string());
}

std::shared_ptr<const System> LoadSystem(OutputDirectory const& directory, std::string const& digest_hex)
{
    auto const document = directory.GetDefinition(digest_hex);
    auto rebuilt = std::make_shared<System>(System::FromCanonicalEncoding(EncodingOfSystemDefinition(document)));
    auto const actual = rebuilt->ContentDigest().Hex();
    if (actual != digest_hex)
        throw std::runtime_error("system definition " + digest_hex
            + ": the rebuilt system's content digest is " + actual
            + "; the document was altered, or its encoding no longer means what it meant when written");
    return InternSystem(rebuilt);
}

} // namespace records
} // namespace bertini
