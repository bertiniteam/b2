//This file is part of Bertini 2.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file eval_helper.hpp

\brief Test-only helper to evaluate a bare function-tree node at a point.

Node-level evaluation has been removed (the straight-line program is the sole evaluator), so tests
that used to set values on nodes and call node->Eval<T>() now evaluate through the SLP.  EvalAt
compiles+runs the expression's program (via EvalExpression) at a point supplied as a name->value
map.  Unlike EvalExpression, it ignores extra entries: callers pass one superset map and each
expression takes only the variables it actually contains, so a single point serves many sub-expressions.
*/

#pragma once

#include <map>
#include <string>
#include <memory>

#include "bertini2/system/eval_expression.hpp"
#include "bertini2/function_tree.hpp"

namespace bertini {
namespace test {

template <typename T>
T EvalAt(std::shared_ptr<node::Node> const& n, std::map<std::string, T> const& known = {})
{
	std::map<std::string, T> point;
	for (auto const& v : node::GatherVariables(n))
		point[v->name()] = known.at(v->name());
	return EvalExpression<T>(n, point);
}

} // namespace test
} // namespace bertini
