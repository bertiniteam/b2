Welcome to PyBertini
====================================

Bertini is software for numerically solving systems of polynomials.  PyBertini is the Python provided for running Bertini.

The main algorithm for numerical algebraic geometry implemented in Bertini is homotopy continuation.  A homotopy is formed, and the solutions to the start system are continued into the solutions for the target system.


.. figure:: images_common/homotopycontinuation_generic_40ppi.png
   :scale: 100 %
   :alt: Homotopy continuation

   Predictor-corrector methods with optional adaptive precision track paths from 1 to 0, solving :math:`f`.




Source
------

PyBertini is distributed with Bertini2, available at `its GitHub repo <https://github.com/bertiniteam/b2>`_.