Welcome to PyBertini
====================================

Bertini is software for numerically solving systems of polynomials.  PyBertini is the Python provided for running Bertini.

Mathematical overview
----------------------

The main algorithm for numerical algebraic geometry implemented in Bertini is homotopy continuation.  A homotopy is formed, and the solutions to the start system are continued into the solutions for the target system.


.. figure:: images_common/homotopycontinuation_generic_40ppi.png
   :scale: 100 %
   :alt: Homotopy continuation

   Predictor-corrector methods with optional adaptive precision track paths from 1 to 0, solving :math:`f`.




Source code
------------

PyBertini is distributed with Bertini2, available at `its GitHub repo <https://github.com/bertiniteam/b2>`_.

The core is written in template-heavy C++, and is exposed to Python through Boost.Python.

Licenses
--------

Bertini2 and its direct components are available under GPL3, with additional clauses in section 7 to protect the Bertini name.  Bertini2 also uses open source softwares, with their own licenses, which may be found in the Bertini2 repository, in the licenses folder.
