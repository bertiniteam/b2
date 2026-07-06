"""Canonical ordering, and opting out -- Bertini 2 tutorial.

Assembled from the ``.. testcode::`` blocks of the "Canonical ordering, and
opting out" tutorial into one runnable program.

Run:  python canonical_ordering.py
"""

import bertini
from bertini import Variable, MonomialOrder


def on_by_default():
    """Commutative reorderings collapse and printed expressions read the conventional way."""
    x, y = Variable('x'), Variable('y')

    assert str(x + y) == str(y + x) == 'x+y'      # commutative: one canonical form (and one node)
    assert str(x * y) == str(y * x) == 'x*y'
    assert str(3 * x**2) == '3*x^2'               # the coefficient leads its monomial
    assert str(x**2 + 2*x - 1) == 'x^2+2*x-1'     # degree-descending; the constant term is last

    assert bertini.canonicalize() is True         # on by default

    return x, y


def choosing_the_monomial_order(x, y):
    """The pluggable, session-global monomial order determines which expressions share a node."""
    bertini.monomial_order(MonomialOrder.Lex)
    lex = str(x**2 + y**3)                         # lexicographic: x before y
    bertini.monomial_order(MonomialOrder.GrevLex)
    grev = str(x**2 + y**3)                        # graded: the higher-degree y^3 leads

    assert lex == 'x^2+y^3'
    assert grev == 'y^3+x^2'

    bertini.monomial_order(MonomialOrder.GrevLex)  # restore the default


def opting_out(x, y):
    """Turn canonicalization off to preserve a deliberately authored operand order."""
    bertini.canonicalize(False)                    # opt out
    try:
        authored = y + x
        assert str(authored) == 'y+x'              # exactly as written, not reordered to x+y
    finally:
        bertini.canonicalize(True)                 # back to the default for everything else


def main():
    x, y = on_by_default()
    choosing_the_monomial_order(x, y)
    opting_out(x, y)


if __name__ == '__main__':
    main()
