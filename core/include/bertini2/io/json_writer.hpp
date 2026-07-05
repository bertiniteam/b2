//This file is part of Bertini 2.
//
//json_writer.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//json_writer.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with json_writer.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file json_writer.hpp

\brief Render a System's structure as JSON: named fields for its parts.

This is the *presentation* view of a system for the structured output directory --
variable groups, path variable, functions, patches -- machine-parseable with one
`json.load`.  Classic (Bertini 1) syntax is deliberately NOT used here: classic is an
input/compatibility format for replicating results in that software, and it cannot
express bertini2's block structure.  The exact, identity-bearing form of a system is
its canonical encoding (`b2sysenc`, ADR-0042); this view is derived and regenerable.
*/

#pragma once

#include <boost/json.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/function_tree/find.hpp"

#include <sstream>
#include <string>

namespace bertini{
	namespace io{

		/**
		\brief A System's structure as a JSON object with fields for its parts.

		Fields: `variable_groups` / `hom_variable_groups` (arrays of arrays of names,
		declaration order), `homogenizing_variables` (names; present in homogenized
		systems), `path_variable` (name, or null), `named_subexpressions`
		(`{name: expression}`, discovered inside the function trees), `functions`
		(expression strings, printed in user coordinates), and `is_patched` (with
		`num_patches` when patched; exact patch coefficients live in the canonical
		encoding, not here).

		\param sys The system to render.
		\return The parts as a boost::json object (presentation only, never identity).
		*/
		inline boost::json::object SystemPartsJson(System const& sys)
		{
			namespace json = boost::json;
			json::object parts;

			auto groups_as_json = [](auto const& groups) {
				json::array out;
				for (auto const& grp : groups)
				{
					json::array names;
					for (auto const& v : grp)
						names.push_back(json::value(v->name()));
					out.push_back(names);
				}
				return out;
			};
			parts["variable_groups"] = groups_as_json(sys.VariableGroups());
			parts["hom_variable_groups"] = groups_as_json(sys.HomVariableGroups());

			{
				json::array homogenizers;
				for (auto const& v : sys.HomogenizingVariables())
					homogenizers.push_back(json::value(v->name()));
				parts["homogenizing_variables"] = homogenizers;
			}

			if (sys.HavePathVariable())
				parts["path_variable"] = sys.GetPathVariable()->name();
			else
				parts["path_variable"] = nullptr;

			auto const functions = sys.NaturalFunctionsAsNodes();
			{
				// named subexpressions are discovered inside the trees (nested included)
				json::object named;
				std::vector<std::shared_ptr<const node::Node>> roots(functions.begin(),
				                                                     functions.end());
				for (auto const& ne : node::Find<node::NamedExpression>(roots))
				{
					std::ostringstream expr;
					expr << *(ne->EntryNode());
					named[ne->name()] = expr.str();
				}
				if (!named.empty())
					parts["named_subexpressions"] = named;
			}
			{
				json::array function_texts;
				for (auto const& f : functions)
				{
					std::ostringstream expr;
					expr << *f;
					function_texts.push_back(json::value(expr.str()));
				}
				parts["functions"] = function_texts;
			}

			parts["is_patched"] = sys.IsPatched();
			if (sys.IsPatched())
				parts["num_patches"] = static_cast<std::int64_t>(sys.NumPatches());

			return parts;
		}

	} // namespace io
} // namespace bertini
