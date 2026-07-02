# This file is part of Bertini 2.
#
# python/bertini/sympy_bridge.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/sympy_bridge.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/sympy_bridge.py.  If not, see <http://www.gnu.org/licenses/>.
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

"""Exact two-way conversion between sympy expressions and bertini function trees.

sympy is an optional dependency of bertini; this module raises a helpful
ImportError if it is missing.

The conversions are EXACT for symbols and exact number types: sympy Rationals
become bertini Rationals (never a float64 round-trip), Integers stay integers
of any size, and Floats travel as strings carrying their precision.  This is
why ``sympy.lambdify`` is deliberately not used -- it lowers ``Rational(1,3)``
to a float literal.

Forward (define in sympy, solve with bertini)::

    import sympy
    from bertini.sympy_bridge import system_from_sympy
    from bertini.nag_algorithm import ZeroDimCauchyAdaptivePrecisionTotalDegree

    x, y = sympy.symbols('x y')
    sys = system_from_sympy([x**2 + y**2 - 1, x + y], [x, y])
    solver = ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)
    solver.solve()
    solver.all_solutions()   # in your variables, comparable with sympy.solve output

Reverse (inspect bertini trees symbolically)::

    from bertini.sympy_bridge import to_sympy
    expr = to_sympy(sys.function(0).root())
"""

try:
    import sympy as _sp
except ImportError as _e:
    raise ImportError(
        "bertini.sympy_bridge requires sympy, which is an optional dependency "
        "of bertini.  install it with: pip install sympy") from _e

import bertini as _pb
import bertini.function_tree as _ft
from bertini.function_tree import operator as _op
from bertini.function_tree import symbol as _sym
from bertini import multiprec as _mp


__all__ = ['from_sympy', 'to_sympy', 'system_from_sympy']


_SYMPY_FUNCTION_MAP = {
    _sp.sin: _ft.sin, _sp.cos: _ft.cos, _sp.tan: _ft.tan,
    _sp.asin: _ft.asin, _sp.acos: _ft.acos, _sp.atan: _ft.atan,
    _sp.exp: _ft.exp, _sp.log: _ft.log,
}

_HYPERBOLICS = (_sp.sinh, _sp.cosh, _sp.tanh, _sp.asinh, _sp.acosh, _sp.atanh)


def _normalize_variables(variables):
    """name -> bertini Variable dict, from a dict or iterable of Variables."""
    if variables is None:
        return {}
    if isinstance(variables, dict):
        return dict(variables)
    return {v.name: v for v in variables}


def from_sympy(expr, variables=None):
    """Convert a sympy expression to a bertini function-tree node, exactly.

    ``variables`` may be an iterable of bertini Variables or a name->Variable
    dict; sympy Symbols are matched to them by name.  Symbols not supplied are
    created fresh, and repeated occurrences share one Variable.

    Raises NotImplementedError for sympy heads bertini has no node for
    (e.g. the hyperbolics -- rewrite first: ``expr.rewrite(sympy.exp)``).
    """
    var_map = _normalize_variables(variables)

    def walk(e):
        if isinstance(e, _sp.Symbol):
            name = e.name
            if name not in var_map:
                var_map[name] = _pb.Variable(name)
            return var_map[name]
        if isinstance(e, _sp.Integer):
            return _sym.Integer(str(e))
        if isinstance(e, _sp.Rational):  # after Integer: Integer is a Rational in sympy
            return _sym.Rational(f'{e.p}/{e.q}')
        if isinstance(e, _sp.Float):
            return _sym.Complex(str(e))
        if e is _sp.pi:
            return _sym.Pi()
        if e is _sp.E:
            return _sym.E()
        if e is _sp.I:
            return _sym.Complex(0, 1)          # the imaginary unit (make_i was Complex(0,1))
        if isinstance(e, _sp.Add):
            result = walk(e.args[0])
            for a in e.args[1:]:
                result = result + walk(a)
            return result
        if isinstance(e, _sp.Mul):
            result = walk(e.args[0])
            for a in e.args[1:]:
                result = result * walk(a)
            return result
        if isinstance(e, _sp.Pow):
            base, exponent = e.args
            if isinstance(exponent, _sp.Integer):
                return walk(base) ** int(exponent)
            return walk(base) ** walk(exponent)
        if isinstance(e, _HYPERBOLICS):
            raise NotImplementedError(
                f"bertini has no node for {type(e).__name__}; rewrite the "
                f"expression first, e.g. expr.rewrite(sympy.exp)")
        # note: sympy has no sqrt head -- sympy.sqrt(x) IS Pow(x, 1/2), handled above
        for head, fn in _SYMPY_FUNCTION_MAP.items():
            if isinstance(e, head):
                return fn(walk(e.args[0]))
        raise NotImplementedError(
            f"no bertini equivalent for sympy node of type {type(e).__name__}: {e}")

    return walk(_sp.sympify(expr))


def system_from_sympy(exprs, variables):
    """Build a bertini System from sympy expressions, with one affine variable group.

    ``variables`` fixes the variable order (iterable of sympy Symbols, bertini
    Variables, or names).  Returns the System; its ``variable_ordering()``
    matches the order given.
    """
    bertini_vars = []
    for v in variables:
        if isinstance(v, _sp.Symbol):
            bertini_vars.append(_pb.Variable(v.name))
        elif isinstance(v, str):
            bertini_vars.append(_pb.Variable(v))
        else:
            bertini_vars.append(v)  # assume already a bertini Variable

    var_map = {v.name: v for v in bertini_vars}
    sys = _pb.System()
    for e in exprs:
        sys.add_function(from_sympy(e, var_map))
    sys.add_variable_group(_pb.VariableGroup(bertini_vars))
    return sys


def _exact_rational(r):
    """bertini multiprec Rational (printed as p/q or integer) -> sympy Rational."""
    return _sp.Rational(str(r))


def to_sympy(node):
    """Convert a bertini function-tree node to a sympy expression, exactly.

    Walks the tree via the introspection API.  Variables map to sympy Symbols
    by name; Integer/Rational leaves convert exactly; Float leaves carry their
    mpfr precision into sympy's.

    Raises NotImplementedError for Differential leaves (trees built by the
    no-argument Jacobian form of differentiate; differentiate with an explicit
    variable instead) and for node types with no sympy equivalent.
    """
    n = node

    if isinstance(n, _op.Sum):
        terms = []
        for i in range(n.num_operands()):
            t = to_sympy(n.operand(i))
            terms.append(t if n.sign(i) else -t)
        return _sp.Add(*terms)
    if isinstance(n, _op.Mult):
        factors = []
        for i in range(n.num_operands()):
            f = to_sympy(n.operand(i))
            factors.append(f if n.mult_or_div(i) else _sp.Pow(f, -1))
        return _sp.Mul(*factors)
    if isinstance(n, _op.Power):
        return _sp.Pow(to_sympy(n.get_base()), to_sympy(n.get_exponent()))
    if isinstance(n, _op.IntegerPower):
        return _sp.Pow(to_sympy(n.operand()), int(n.exponent))
    if isinstance(n, _op.Negate):
        return -to_sympy(n.operand())
    if isinstance(n, _op.Sqrt):
        return _sp.sqrt(to_sympy(n.operand()))
    if isinstance(n, _op.Exp):
        return _sp.exp(to_sympy(n.operand()))
    if isinstance(n, _op.Log):
        return _sp.log(to_sympy(n.operand()))
    if isinstance(n, _op.Sin):
        return _sp.sin(to_sympy(n.operand()))
    if isinstance(n, _op.Cos):
        return _sp.cos(to_sympy(n.operand()))
    if isinstance(n, _op.Tan):
        return _sp.tan(to_sympy(n.operand()))
    if isinstance(n, _op.ArcSin):
        return _sp.asin(to_sympy(n.operand()))
    if isinstance(n, _op.ArcCos):
        return _sp.acos(to_sympy(n.operand()))
    if isinstance(n, _op.ArcTan):
        return _sp.atan(to_sympy(n.operand()))

    if isinstance(n, _sym.Differential):
        raise NotImplementedError(
            "this tree contains Differential leaves (built by the no-argument, "
            "Jacobian form of differentiate); differentiate with an explicit "
            "variable instead, e.g. f.differentiate(x)")
    if isinstance(n, _sym.Variable):
        return _sp.Symbol(n.name)
    if isinstance(n, _sym.Pi):
        return _sp.pi
    if isinstance(n, _sym.E):
        return _sp.E
    if isinstance(n, _sym.Integer):
        return _sp.Integer(int(str(n.value())))
    if isinstance(n, _sym.Rational):
        re = _exact_rational(n.value_real())
        im = _exact_rational(n.value_imag())
        return re if im == 0 else re + _sp.I * im
    if isinstance(n, _sym.Complex):
        v = n.value()
        digits = v.real.precision
        # repr is full-precision (str truncates to ostream's default 6 digits)
        re = _sp.Float(repr(v.real), digits)
        if v.imag == _mp.real_mp(0):
            return re
        return re + _sp.I * _sp.Float(repr(v.imag), digits)

    raise NotImplementedError(
        f"no sympy equivalent for bertini node of type {type(n).__name__}: {n}")
