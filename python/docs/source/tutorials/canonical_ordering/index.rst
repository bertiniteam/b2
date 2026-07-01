🧮 Canonical ordering, and opting out
*************************************************

.. testsetup:: *

   import bertini
   from bertini import Variable, MonomialOrder

Bertini hash-conses the function tree: identical subexpressions become a single shared node.
For that sharing to fire on expressions you build *separately*, the library also **canonically
orders** the operands of every sum and product -- so ``x + y`` and ``y + x`` are not just equal,
they are the **same node**. This tutorial shows what that ordering does, how to choose it, and how
to turn it off when you want to keep the exact operand order you wrote.

On by default
=============

Commutative reorderings collapse, coefficients sort to the front of a monomial, and the constant
term of a sum stays last -- so printed expressions read the conventional way.

.. testcode::

   x, y = Variable('x'), Variable('y')

   assert str(x + y) == str(y + x) == 'x+y'      # commutative: one canonical form (and one node)
   assert str(x * y) == str(y * x) == 'x*y'
   assert str(3 * x**2) == '3*x^2'               # the coefficient leads its monomial
   assert str(x**2 + 2*x - 1) == 'x^2+2*x-1'     # degree-descending; the constant term is last

You can ask whether canonicalization is on, and toggle it, with :func:`bertini.canonicalize`.

.. testcode::

   assert bertini.canonicalize() is True         # on by default

Choosing the monomial order
===========================

The order is a pluggable monomial order -- ``Lex``, ``RevLex``, or ``GrevLex`` (the default,
graded so higher total degree leads). It is a session-global setting, because it determines which
expressions share a node.

.. testcode::

   bertini.monomial_order(MonomialOrder.Lex)
   lex = str(x**2 + y**3)                         # lexicographic: x before y
   bertini.monomial_order(MonomialOrder.GrevLex)
   grev = str(x**2 + y**3)                        # graded: the higher-degree y^3 leads

   assert lex == 'x^2+y^3'
   assert grev == 'y^3+x^2'

   bertini.monomial_order(MonomialOrder.GrevLex)  # restore the default

Opting out
==========

Canonicalization is the right default -- it maximizes sharing and keeps the compiled
straight-line program small. But sometimes you have arranged the operands of an expression
**deliberately**, and you want the library to leave that arrangement alone. Turn canonicalization
off while you build that expression, and the authored operand order is preserved:

.. testcode::

   bertini.canonicalize(False)                    # opt out
   try:
       authored = y + x
       assert str(authored) == 'y+x'              # exactly as written, not reordered to x+y
   finally:
       bertini.canonicalize(True)                 # back to the default for everything else

Why would you do this? Floating-point addition and multiplication are not associative: the order
in which operands are combined changes the rounding. If you have written an expression so that a
catastrophically-cancelling pair is combined first (or so that a product avoids an intermediate
overflow), canonical reordering could undo that arrangement. Opting out hands you back exact
control of the operand order for the numerically-sensitive part of a model, while everything you
do *not* opt out of stays canonicalized -- and keeps sharing.

.. note::

   The toggle is session-global, so the idiom is "opt out, build the sensitive expression, opt
   back in" (as above). Everything built while it is off keeps its authored order and is a
   distinct node from the canonical form; everything built while it is on is canonicalized and
   shared.
