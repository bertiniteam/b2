//This file is part of Bertini 2.
//
//reintern.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//reintern.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with reintern.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file reintern.hpp

\brief Re-intern deserialized node DAGs into the live hash-consing tables (ADR-0042).

Deserialization deliberately constructs nodes WITHOUT going through Make/Intern
(node.hpp), so a freshly loaded tree does not share with structurally equal live nodes
-- and on a System-intern miss the loaded nodes would fork the intern universe: nothing
built afterward, and nothing loaded from a second archive, would unify with them.
Reintern() repairs that: a memoized post-order rebuild of the loaded DAG through the
ordinary interning factories, so every subtree collapses onto its live canonical
representative (or becomes it).

Runs ONLY at deserialization time.  The pointer-keyed memo preserves the archive's
internal sharing (each loaded object rebuilt once; the pass is linear).  Rebuilding
re-runs construction-time canonicalization, so it assumes the session canonicalization
settings match those the tree was authored under -- which the System digest header
records (ADR-0042).  Rebuilt nodes acquire fresh memoized state (structural hash,
compiled-evaluator cache).  Content is never changed: reinterning a tree leaves its
canonical encoding, and hence any digest over it, bit-identical.
*/

#pragma once

#include <map>
#include <memory>

#include "bertini2/function_tree/node.hpp"

namespace bertini {
namespace node {

/// Memo for one re-interning pass: loaded node -> its interned replacement.  Share one
/// memo across all roots of a loaded object (e.g. a whole System), so cross-root sharing
/// in the archive is preserved in the rebuilt DAG.
using ReinternMemo = std::map<Node const*, std::shared_ptr<Node>>;

/**
\brief Rebuild the DAG rooted at `root` through the interning factories (memoized,
post-order), returning the canonical interned equivalent.

\param root The (typically just-deserialized) root to re-intern.
\param memo The pass's memo; share it across the roots of one loaded object.
\return The interned equivalent (a live pre-existing node when one matches).

Throws on a node kind unknown to the rebuilder (fail loud -- an unrebuildable node
must never silently remain un-interned).
*/
std::shared_ptr<Node> Reintern(std::shared_ptr<Node> const& root, ReinternMemo& memo);

/**
\brief Convenience overload with a private memo (single-root use).

\param root The root to re-intern.
\return The interned equivalent.
*/
std::shared_ptr<Node> Reintern(std::shared_ptr<Node> const& root);

} // namespace node
} // namespace bertini
