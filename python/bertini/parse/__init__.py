# This file is part of Bertini 2.
# 
# python/bertini/parse/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/parse/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/parse/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  silviana amethyst
#  UWEC
#  Spring 2018
# 





"""
Parsing functions, taking strings and producing various other things
"""

import re as _re

from bertini._pybertini import parse as _pybparse
from bertini._pybertini.parse import *

# The native SystemParser wants JUST the input body (variable groups + functions); it
# rejects a leading CONFIG...END; block and does not want the INPUT/END; separators that
# `System.to_classic_input()` itself emits -- so `parse.system(a_system.to_classic_input())`
# fails and callers were forced to hand-split the text.  Wrap it to tolerate both, so the
# round-trip just works and no caller has to know the parser's quirks.

_native_system = _pybparse.system


def _input_body(text):
    """The parseable body of a classic Bertini input: drop a leading CONFIG...END; block
    and, if the remainder is wrapped in INPUT...END;, return just the inside.  A raw body
    (no CONFIG, no INPUT wrapper) is returned unchanged, so this never breaks input the
    native parser already accepted."""
    s = text
    m = _re.search(r'\bCONFIG\b.*?\bEND\s*;', s, flags=_re.IGNORECASE | _re.DOTALL)
    if m:
        s = s[:m.start()] + s[m.end():]
    m = _re.search(r'\bINPUT\b(.*?)\bEND\s*;', s, flags=_re.IGNORECASE | _re.DOTALL)
    if m:
        return m.group(1).strip()
    return s.strip()


def system(text):
    """Parse a classic Bertini system from ``text``, tolerating a leading CONFIG section
    and INPUT/END; separators (so ``parse.system(sys.to_classic_input())`` round-trips).
    """
    return _native_system(_input_body(text))


__all__ = dir(_pybparse)

