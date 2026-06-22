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

The work is delegated to the single evaluation engine (the compiled SLP).  The expression's
variables are discovered and ordered by name; a one-function program is compiled for the
expression and **memoized on the (immutable, hash-consed) node**, so repeated evaluation of the
same expression reuses it rather than recompiling.  Each call runs the shared program against
its own working memory at the requested point.
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
	by name --- the session-global canonical order), a one-function program is compiled for the
	expression (memoized on the node) and run at the point.  Variables are matched to values by
	name, which is unambiguous because variables are canonical by name across the whole session.

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

		// Compile the expression's evaluator once and memoize it on the (immutable, hash-consed)
		// node, so repeated evaluation of the same expression reuses the compiled program rather
		// than rebuilding a throwaway System and recompiling an SLP every call (ADR-0027).
		// The compiled program holds no node pointers (its constants are value recipes), so caching
		// it on the node creates no reference cycle.
		auto cached = std::static_pointer_cast<const StraightLineProgram>(expr->EvalProgram());
		if (!cached)
		{
			System sys;
			sys.AddFunction(expr);
			sys.AddVariableGroup(vars);
			cached = std::make_shared<const StraightLineProgram>(sys);
			expr->SetEvalProgram(cached);
		}

		// Per-call working copy: it shares the immutable compiled Program but gets its own
		// evaluation Memory, so concurrent evaluations of the same cached expression do not race.
		StraightLineProgram slp = *cached;

		if constexpr (!std::is_same<T, dbl>::value)
			if (vars.size() > 0)
				slp.precision(Precision(point));

		slp.Eval(point);
		return slp.template GetFuncVals<T>()(0);
	}

} // namespace bertini
