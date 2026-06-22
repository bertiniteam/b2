//This file is part of Bertini 2.
//
//include/bertini2/system/eval_expression.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/system/eval_expression.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/system/eval_expression.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


/**
\file include/bertini2/system/eval_expression.hpp

\brief Evaluate a bare expression tree at a point, without the caller owning a System.

The work is delegated to the single evaluation engine (the compiled SLP): a throwaway
single-function adapter System is built around the expression, its variables discovered
and ordered by name, and its SLP is run at the requested point.  No new evaluation path is
introduced --- this is a convenience wrapper over `System::Eval`.
*/

#pragma once

#include <map>
#include <string>
#include <stdexcept>
#include <type_traits>

#include "bertini2/system/system.hpp"
#include "bertini2/function_tree/gather.hpp"

namespace bertini {

	/**
	\brief Evaluate an expression tree at a point supplied as a variable-name -> value map.

	The expression's variables are discovered with `node::GatherVariables` (distinct, sorted
	by name --- the session-global canonical order), wrapped in a one-function adapter System,
	and evaluated through that System's compiled SLP.  Variables are matched to values by name,
	which is unambiguous because variables are canonical by name across the whole session.

	Every variable appearing in the expression must have a supplied value, and every supplied
	name must appear in the expression; either kind of mismatch throws.  The latter is a
	guard against typos --- a value silently going nowhere is almost always a mistake.

	\tparam T The evaluation number type (`dbl` or `mpfr_complex`).
	\param expr The expression to evaluate.
	\param variable_values A map from variable name to the value to substitute.
	\return The value of the expression at the given point.
	*/
	template<typename T>
	T EvalExpression(std::shared_ptr<node::Node> const& expr,
	                 std::map<std::string, T> const& variable_values)
	{
		auto vars = node::GatherVariables(expr); // distinct, sorted by name

		// Guard against typos: every supplied name must be a variable of the expression.
		for (auto const& named_value : variable_values)
		{
			bool appears = false;
			for (auto const& v : vars)
				if (v->name() == named_value.first) { appears = true; break; }
			if (!appears)
				throw std::runtime_error("value supplied for '" + named_value.first +
					"', which is not a variable of the expression being evaluated");
		}

		// Bind values into the System's variable order (which is the name order of `vars`).
		Vec<T> point(static_cast<Eigen::Index>(vars.size()));
		for (size_t i = 0; i < vars.size(); ++i)
		{
			auto it = variable_values.find(vars[i]->name());
			if (it == variable_values.end())
				throw std::runtime_error("no value supplied for variable '" + vars[i]->name() +
					"' when evaluating an expression");
			point(static_cast<Eigen::Index>(i)) = it->second;
		}

		System sys;
		sys.AddFunction(expr, "expr");
		sys.AddVariableGroup(vars);

		if constexpr (!std::is_same<T, dbl>::value)
			if (vars.size() > 0)
				sys.precision(Precision(point));

		return sys.template Eval<T>(point)(0);
	}

} // namespace bertini
