//This file is part of Bertini 2.
//
//canonical.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical.hpp

\brief Canonical ordering of commutative-operator operands (reorder-only).

When canonicalization is enabled, the operands of a Sum / Mult are sorted into a
deterministic, content-based canonical order BEFORE the node is interned, so that
structurally-equal-up-to-reordering expressions (x+y and y+x) collapse to one shared node.

The order is a pluggable monomial order -- Lex / RevLex / GrevLex -- on the operands'
exponent vectors (via Node::MultiDegree), built on the global variable order BY NAME
(node::GatherVariables already returns variables sorted alphabetically).  Non-polynomial
operands sort after the polynomial ones, with a printed-form tie-break to stay total.
This sorts operands but never combines them (no like-term folding; a full polynomial normal
form is future work).

Nothing here changes Node::Hash()/IsSame(): canonicalization is a normalization applied at
construction, and the existing order-sensitive predicates then see the normalized order.
*/

#pragma once

#include <vector>
#include <memory>

#include "bertini2/function_tree/node.hpp"

namespace bertini {
namespace node {

/// Pluggable monomial orders for canonicalization.  The default is deferred; GrevLex (the
/// usual computational-algebra default, and graded => degree-descending, matching how
/// polynomials are conventionally written) is the placeholder.
enum class MonomialOrder { Lex, RevLex, GrevLex };

/// The monomial order currently used for canonicalization (session-global -- it must be
/// consistent session-wide, since it determines interned identity).
MonomialOrder CurrentMonomialOrder();
void SetMonomialOrder(MonomialOrder order);

/// Whether Sum/Mult construction canonicalizes operand order by default.  (Per-expression
/// opt-out is a future refinement; for now this is the session switch.)
bool CanonicalizeByDefault();
void SetCanonicalizeByDefault(bool on);

/**
\brief Sort an n-ary operator's (operand, flag) pairs into canonical order, in place.

\param operands the operator's operands (mutated into canonical order)
\param flags     the parallel signs (Sum) or multiply/divide flags (Mult), reordered in tandem
\param multiplicands_first  for a Mult, keep flag==true (multiplicands) ahead of flag==false
       (divisors), so a divisor never ends up first (the eval path expects a leading
       multiplicand).  Pass false for a Sum (the +/- sign does not affect term order).

No-op when canonicalization is disabled or there are fewer than two operands.  Called from
the Sum/Mult constructors, so the node is canonical before Make() interns it.
*/
void CanonicalizeNaryOperands(std::vector<std::shared_ptr<Node>>& operands,
                              std::vector<bool>& flags,
                              bool multiplicands_first);

} // namespace node
} // namespace bertini
