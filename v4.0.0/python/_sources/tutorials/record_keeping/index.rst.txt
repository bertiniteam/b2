The record-keeping system
***************************

Solving systems using numerical algebraic geometry can be expensive.  So we implemented a recording system that's on
by default, so when you solve a system, the results are automatically saved to disk
in a way that lets re-solves of the exact same system with the same settings
recall those solutions.  So solving writes a durable, plain-text structured output directory, and
recalls from it so reruns are instant and crashes resume about where they left off.
These tutorials teach about the recording system:
the automatic records every solve keeps, and how chained solves thread provenance through
them.  I hope you like it!

.. toctree::
   :maxdepth: 1

   automatic_record_keeping/index
   chained_homotopies/index
