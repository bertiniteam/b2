//This file is part of Bertini 2.
//
//print_dialect.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//print_dialect.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with print_dialect.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file print_dialect.hpp

\brief The dialect a function tree is printed in -- who the text is for.

Printing is a family, one member per audience (ADR-0059):

- **Classic**: Bertini 1 input.  Powers are `^`; a complex constant is `(re+im*I)`; every
  constant carries all of its digits.  This is what `System::to_classic_input` and the CLI
  write, and what the classic parser reads back.
- **PythonReadable**: what `str(node)` shows in Python.  Powers are `**`; constants keep the
  same shapes as Classic (a bare real, `(re+im*I)`); readable first, exact where that costs
  nothing.
- **PythonExact**: what `repr(node)` shows.  `eval` of it in the `bertini` namespace rebuilds
  an equal node at full precision: a multiprecision real is `real_mp('digits', precision)`, a
  multiprecision complex constant is `Complex('re', 'im', precision)`, a non-integer rational
  is `Rational('p/q')`.  Integers, and values a Python literal represents exactly, stay bare.

The dialect rides on the output stream (an `std::ios_base` slot), so the recursive node
printers need no extra parameter and `operator<<` keeps working; a stream that nobody has
set a dialect on prints Classic, which is what every pre-existing caller expected.
Canonicalization (the content digest) uses none of these: identity is decided on the
canonical encoding, never on printed text.
*/

#pragma once

#include <ios>

namespace bertini {
namespace node {

/// \brief Who a printed expression is for; selects the operator and constant spellings.
enum class PrintDialect
{
    Classic = 0,          ///< Bertini 1 input: `^`, `(re+im*I)`, all digits.
    PythonReadable = 1,   ///< Python `str`: `**`, the same constant shapes, readable first.
    PythonExact = 2       ///< Python `repr`: `**`, exact constant spellings that `eval` rebuilds at full precision.
};

/// \brief The stream-state slot the dialect is kept in (allocated once per process).
inline int PrintDialectSlot()
{
    static int const slot = std::ios_base::xalloc();
    return slot;
}

/// \brief The dialect a stream prints in; Classic unless SetDialect was called on it.
inline PrintDialect DialectOf(std::ios_base& stream)
{
    return static_cast<PrintDialect>(stream.iword(PrintDialectSlot()));
}

/// \brief Choose the dialect everything subsequently printed to `stream` is written in.
inline void SetDialect(std::ios_base& stream, PrintDialect dialect)
{
    stream.iword(PrintDialectSlot()) = static_cast<long>(dialect);
}

/// \brief The power operator in a dialect: `^` for Bertini 1, `**` for Python.
inline char const* PowerSymbol(PrintDialect dialect)
{
    return dialect == PrintDialect::Classic ? "^" : "**";
}

} // namespace node
} // namespace bertini
