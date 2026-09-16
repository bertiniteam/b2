//This file is part of Bertini 2.
//
//canonical_decoding.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_decoding.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_decoding.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical_decoding.hpp

\brief The reader for the canonical exact text encoding of expression trees -- the
inverse of canonical_encoding.hpp (ADR-0042).

The canonical encoding is the digest preimage the records archive stores for every
system and homotopy.  This reader turns that text back into an interned node DAG, so an
archived system can be rebuilt with its content digest as the proof of fidelity: the
decoded object encodes to the very same text.  It is a plain recursive-descent reader
over the grammar sketched in canonical_encoding.hpp; every node is rebuilt through its
ordinary Make factory, so the result lands in (or comes from) the live intern tables.

Back-references (`#<index>`) resolve through the same first-encounter numbering the
encoder used: a node receives its index when its opening parenthesis is read, before
its children, exactly as the encoder numbers a node before descending.
*/

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "bertini2/function_tree/node.hpp"

namespace bertini {
namespace node {

/**
\brief A read position over one canonical-encoding text, with the small lexical
helpers the grammar needs (literals, words, netstrings, numbers).

The cursor does not own the text: keep the string alive for the cursor's lifetime.
Every helper throws std::runtime_error with the byte offset and a snippet of the text
when the input does not match, so a malformed or truncated encoding fails loudly at
the first wrong byte instead of producing a wrong object.
*/
class DecodingCursor
{
public:
	/// \brief Start reading `text` from its first byte.
	explicit DecodingCursor(std::string const& text) : text_(text), pos_(0) {}

	/// \brief Whether every byte has been consumed.
	bool AtEnd() const { return pos_ >= text_.size(); }

	/// \brief The byte offset of the next unread character.
	std::size_t Offset() const { return pos_; }

	/// \brief The next unread character (throws at end of text).
	char Peek() const;

	/// \brief Consume exactly the character `c`, or throw.
	void Expect(char c);

	/// \brief Consume exactly the literal `lit`, or throw.
	void Expect(std::string_view lit);

	/// \brief Consume `c` if it is next; report whether it was.
	bool TryConsume(char c);

	/// \brief Skip any run of spaces and newlines.
	void SkipWhitespace();

	/// \brief Read a maximal run of characters that are not whitespace or a parenthesis.
	std::string ReadWord();

	/// \brief Read a netstring name: `<byte-length>:<bytes>`.
	std::string ReadNetstring();

	/// \brief Read a non-negative decimal integer.
	unsigned long long ReadUnsigned();

	/// \brief Read a decimal integer with an optional sign.
	long long ReadInt();

	/// \brief Throw std::runtime_error describing `what` at the current position.
	[[noreturn]] void Fail(std::string const& what) const;

private:
	std::string const& text_;   ///< The encoding being read (not owned).
	std::size_t pos_;           ///< Offset of the next unread byte.
};

/**
\brief The back-reference table of one decoding pass: the nodes in first-encounter
order, so `#<index>` resolves to the node the encoder numbered `index`.

Share one context across every root of an encoding unit (a whole System, with the
operand systems nested inside its blocks), mirroring the EncodingContext the encoder
shared across those same roots.
*/
struct DecodingContext
{
	/// The node numbered by each index, in encounter order.
	std::vector<std::shared_ptr<Node>> by_index;
};

/**
\brief Decode one node (and everything under it) from the cursor's position.

\param cursor The read position; on return it sits just past the node's closing
       parenthesis (or past a back-reference).
\param ctx The back-reference table shared across the encoding unit.
\return The rebuilt, interned node.

Throws std::runtime_error on any byte that does not fit the grammar, on an unknown node
kind, and on a back-reference to an index not yet encountered.
*/
std::shared_ptr<Node> DecodeCanonical(DecodingCursor& cursor, DecodingContext& ctx);

/**
\brief Decode a whole tree from its canonical encoding text (fresh context).

The entire text must be one node: trailing bytes are an error.

\param text The canonical encoding, as CanonicalEncoding(root) produced it.
\return The rebuilt, interned node.
*/
std::shared_ptr<Node> DecodeCanonicalTree(std::string const& text);

} // namespace node
} // namespace bertini
