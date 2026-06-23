//This file is part of Bertini 2.
//
//Bertini 2 is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Bertini 2 is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with Bertini 2. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2026 by Bertini2 developers
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file classic_writer.hpp

\brief Emit a System (and, with a config, a whole solve) as a Bertini 1 *classic* input file.

The classic parser (io/parsing) reads Bertini 1 input files; this is its inverse.  Emitting a
`CONFIG ... END;\nINPUT ... END;` file lets the *same* problem be run in Bertini 1 for
cross-validation -- e.g. to check whether a path-crossing or a root count reproduces across
implementations (modulo each one's random start system).

The emitted text is round-trip-safe: parsing it back with bertini::System's classic constructor
yields an equivalent system.
*/

#pragma once

#include "bertini2/system/system.hpp"
#include "bertini2/function_tree/find.hpp"

#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace bertini{
	namespace classic{

		/**
		\brief Emit the INPUT-section body of a system in classic syntax: the variable groups, any
		named subexpressions (as `name = expr;` subfunctions, defined before the functions that use
		them), and the functions (`function f0,f1,...;` then `f0 = ...;`).

		No `INPUT`/`END;` wrapper -- WriteClassicInput adds that.  Functions are named `f0, f1, ...`
		positionally (Bertini 2 functions carry no user name after parsing).
		*/
		inline void EmitSystem(std::ostream& out, System const& sys)
		{
			auto emit_groups = [&out](char const* keyword, auto const& groups){
				for (auto const& grp : groups)
				{
					out << keyword << " ";
					for (size_t i = 0; i < grp.size(); ++i)
						out << (i ? ", " : "") << *grp[i];
					out << ";\n";
				}
			};
			emit_groups("variable_group", sys.VariableGroups());
			emit_groups("hom_variable_group", sys.HomVariableGroups());

			auto functions = sys.NaturalFunctionsAsNodes();

			// Named subexpressions are not stored separately; they are discovered inside the function
			// trees (nested ones included).  Emit each as `name = expr;` before the functions, the way
			// the classic parser expects a subfunction to be defined ahead of its use.
			{
				std::vector<std::shared_ptr<const node::Node>> roots(functions.begin(), functions.end());
				for (auto const& ne : node::Find<node::NamedExpression>(roots))
					out << ne->name() << " = " << ne->EntryNode() << ";\n";
			}

			out << "function ";
			for (size_t i = 0; i < functions.size(); ++i)
				out << (i ? ", " : "") << "f" << i;
			out << ";\n";
			for (size_t i = 0; i < functions.size(); ++i)
				out << "f" << i << " = " << functions[i] << ";\n";
		}

		/// \brief The INPUT-section body as a string.  \see EmitSystem.
		inline std::string SystemToClassic(System const& sys)
		{
			std::ostringstream ss;
			EmitSystem(ss, sys);
			return ss.str();
		}

	} // namespace classic
} // namespace bertini
