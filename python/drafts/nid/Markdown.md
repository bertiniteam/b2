### Drafts of python versions of code that will make it's way into core

`SliceMover.ipynb` contains a class that can create homotopies and calculate the slice move on a system using the bertini.nag_algorithm.moving_homotopy(). It contains example driver code as well as a sample graph, other then the driver code, and it only requires access to the bertini library. 

`RegenCascade.ipynb` contains the beginning of the Regeneration Cascade to generate witness points for the system. It requires junk removal to be added still, but it'll start generating the points needed for the next step. Currently it generates start points from codimension 1, and then I'm working on codim 2 and beyond.
