👋 Welcome to Bertini 2
====================================

Bertini 2 is software for numerically solving systems of polynomials. 

🗺 Mathematical overview
----------------------------------

The main algorithm for numerical algebraic geometry implemented in Bertini is homotopy continuation.  A homotopy is formed, and the solutions to the start system are continued into the solutions for the target system.


.. figure:: images_common/homotopycontinuation_generic.png
   :scale: 100 %
   :height:  285 px
   :width: 400 px
   :alt: Homotopy continuation

   Predictor-corrector methods with optional adaptive precision track paths from 1 to 0, solving :math:`f`.

The definitive resource for Bertini 1 is the book :cite:`bertinibook`.  While the way we interact with Bertini changes from version 1 to version 2, particularly when using the Python library, the algorithms remain fundamentally the same.  So do most of the ways to change settings for the path trackers, etc.  We believe that embracing the flexibility of Python allows for much greater flexibility.  It also will relieve the user from the burden of input and output file writing and parsing.  Instead, computed results are returned directly to the user.  

Consider checking out the :ref:`🔦 Tutorials <tutorials>`.


⛲️ Source code
---------------------------

The Bertini 2 source code is available at `its GitHub repo <https://github.com/bertiniteam/b2>`_.

The core is written in template-heavy C++, and is exposed to Python through Boost.Python.

⚖️ Licenses
------------------

Bertini2 and its direct components are available under GPL3, with additional clauses in section 7 to protect the Bertini name.  Bertini2 also uses open source softwares, with their own licenses, which may be found in the Bertini2 repository, in the licenses folder.
