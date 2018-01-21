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

The definitive resource for Bertini 1 is the book :cite:`bertinibook`.  While the way we interact with Bertini changes from version 1 to version 2, particularly when using PyBertini, the algorithms remain fundamentally the same.  So do most of the ways to change settings for the path trackers, etc.  We believe that embracing the flexibility of Python3 with PyBertini allows for much greater flexibility.  It also will relieve the user from the burden of input and output file writing and parsing.  Instead, computed results are returned directly to the user.  

Consider checking out the :ref:`tutorials`.


Source code
------------

PyBertini is distributed with Bertini2, available at `its GitHub repo <https://github.com/bertiniteam/b2>`_.

The core is written in template-heavy C++, and is exposed to Python through Boost.Python.

Licenses
--------

Bertini2 and its direct components are available under GPL3, with additional clauses in section 7 to protect the Bertini name.  Bertini2 also uses open source softwares, with their own licenses, which may be found in the Bertini2 repository, in the licenses folder.
