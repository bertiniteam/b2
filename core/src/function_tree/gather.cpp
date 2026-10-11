//This file is part of Bertini 2.
//
//src/function_tree/gather.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/gather.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/gather.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/function_tree/gather.hpp"
#include "bertini2/function_tree.hpp"

#include <algorithm>

namespace bertini {
namespace node {

    // GatherVariables is what Find<Variable> finds -- distinct, sorted by name -- but read from
    // each node's degree record (Node::Variables), which is made once per node from its
    // operands' (#505): canonicalization asks it of every operand of every sum and product.

    VariableGroup GatherVariables(std::shared_ptr<const Node> const& n)
    {
        if (!n)
            return {};
        return n->Variables();
    }

    VariableGroup GatherVariables(std::vector<std::shared_ptr<Node>> const& functions)
    {
        VariableGroup out;
        for (auto const& f : functions)
            if (f)
                out.insert(out.end(), f->Variables().begin(), f->Variables().end());
        std::sort(out.begin(), out.end(),
                  [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b) { return a->name() < b->name(); });
        out.erase(std::unique(out.begin(), out.end(),
                              [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b) { return a.get() == b.get(); }),
                  out.end());
        return out;
    }

} // namespace node
} // namespace bertini
