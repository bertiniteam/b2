//This file is part of Bertini 2.
//
//canonical_encoding.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_encoding.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_encoding.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical_encoding.hpp

\brief The canonical exact text encoding of expression trees (ADR-0042).

This is the identity substrate for persistent content digests: a deterministic,
bit-exact, S-expression-like encoding of a node DAG, stable across runs, compilers,
and library versions.  It deliberately does NOT reuse:

- `Node::Hash()` / typeid -- per-process only;
- the classic printer (ADR-0012) -- a cosmetic, test-pinned surface whose tweaks must
  not rotate persistent digests;
- Boost.Serialization archives -- their bytes embed boost-version-dependent framing.

Format sketch (kind tags are fixed strings, never typeid):

- leaves: `(var 1:x)`, `(diff 1:x)`, `(int -7)`, `(rat 1/3 0)`, `(cplx 30 <re> <im>)`
  (a Complex literal carries its stored precision -- value-equal literals at different
  precisions are distinct identities), `(pi)`, `(e)`;
- operators: `(sum +<child> -<child> ...)` / `(mul *<child> /<child> ...)` in the order
  the session canonicalization produced, `(neg <c>)`, `(pow <base> <exp>)`,
  `(ipow <n> <c>)`, `(sqrt <c>)`, `(exp <c>)`, `(log <c>)`, `(sin <c>)` etc.,
  `(named 1:f <c>)`;
- names are netstrings (`<byte-length>:<utf8-bytes>`), so any name -- emoji included --
  encodes unambiguously;
- back-references: every node gets an index in traversal order at its first encounter;
  a revisit emits `#<index>`.  On interned (hash-consed) DAGs equal subtrees ARE the
  same object, so the encoding is canonical for equal content and stays linear-size
  under heavy sharing.

Any change to this format is a digest-breaking change: bump the system-encoding
version header (`b2sysenc/<n>`, see system content identity) and the golden-digest
fixture in the same commit.
*/

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include "bertini2/function_tree/node.hpp"

namespace bertini {
namespace node {

/**
\brief Traversal state for one canonical-encoding pass: the back-reference table.

One context spans one encoding unit (a whole System, or a single tree), so nodes
shared ACROSS functions back-reference deterministically within that unit.
*/
struct EncodingContext
{
	/// First-encounter index of each node visited so far, keyed by object identity.
	std::map<Node const*, unsigned> seen;
	/// The index the next first-encountered node will receive (traversal order).
	unsigned next_index = 0;
};

/**
\brief Append the canonical encoding of the DAG rooted at `root` to `out`.

\param root The root of the (interned) expression DAG to encode.
\param out The stream the encoding is appended to.
\param ctx The back-reference table; share one context across all roots of an
       encoding unit, so cross-root sharing encodes deterministically.

Throws on a node kind unknown to the encoder (fail loud: an unencodable node must
never silently produce a digest).
*/
void EncodeCanonical(std::shared_ptr<Node> const& root, std::ostream& out, EncodingContext& ctx);

/**
\brief The canonical encoding of a single tree as a string (fresh context).

\param root The root of the (interned) expression DAG to encode.
\return The canonical encoding text.
*/
std::string CanonicalEncoding(std::shared_ptr<Node> const& root);

} // namespace node
} // namespace bertini
