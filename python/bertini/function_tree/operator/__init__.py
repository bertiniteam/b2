# This file is part of Bertini 2.
#
# python/bertini/function_tree/operator/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/function_tree/operator/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/function_tree/operator/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#  silviana amethyst
#  University of Wisconsin - Eau Claire
#  Spring 2026
#

"""Operator node types (Sum, Mult, Power, ...), for inspecting function trees.

Unary, Nary, and Trig are abstract but deliberately kept: they make
isinstance-based tree walking possible (e.g. any unary operator exposes
operand()).
"""

from bertini._pybertini.function_tree import operator as _pybop
from bertini._pybertini.function_tree.operator import *

del AbstractOp

_ABSTRACT = {'AbstractOp'}
__all__ = [n for n in dir(_pybop) if n not in _ABSTRACT]
