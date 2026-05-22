
# Quick links

- [Wiki](https://github.com/bertiniteam/b2/wiki)

---

# Overview

The solution of arbitrary polynomial systems is an area of active research, and has many applications in math, science and engineering.  This software, Bertini 2, is a complete re-implementation of [Bertini 1](https://bertini.nd.edu) from C into C++/Python.

The theoretical basis for the solution of polynomials with Bertini is a theorem which bounds the number of solutions a system may have. It sits together with the numerical computational tool of "homotopy continuation". the act of "continuing" from one system into another through a "homotopy", as depicted in the below diagram:

<img src="https://raw.githubusercontent.com/bertiniteam/b2/develop/doc_resources/images/homotopycontinuation_generic.png" alt="homotopy continuation"
	title="homotopy continuation" width="400" height="285"/>

---

# Current capabilites

Bertini 2 currently has implemented the foundations of Numerical Algebraic Geometry.  Development is ongoing, and here's what we have so far:

- C++ functions and types, with Python bindings.
- Through Python, runtime scriptable construction of systems and interactivity with their zero-dimensional solutions.
- Construction of multivariate polynomial and non-polynomial systems.
- Evaluation of systems and their Jacobians in double and arbitrary multiple precision, using two different methods.
- Construction of the Total Degree and Multihomogeneous start systems.
- Construction of homotopies (they're just systems with path variables defined).
- Tracking of a start point x_0, corresponding to a particular time $t_0 \in \mathbb{C}^n$ in a homotopy $H$, from $t_0$ to $t_1$.
- Running of the Power Series and Cauchy endgames, in double, multiple, and adaptive precision.

Development is ongoing, and we want your help!

---

# Missing functionality

* Parallel solving
* Numerical irreducible decomposition
* Membership testing
* and other algorithms

Users wanting a more developed implementation are recommended to use [Bertini 1](https://bertini.nd.edu) or [homotopycontinuation.jl](https://www.juliahomotopycontinuation.org/), or one of the other packages implementing the theory.  

---

# Installation

## Pre-built wheels -- the way to go!

The Python package `bertini2` provides pre-built wheels for Linux, macOS, and Windows. 


```bash
pip install bertini2
```

Once it's installed, you `import bertini`

* Linux: Python 3.9-3.14
* MacOS (Apple Silicon): Python 3.9-3.14
* MacOS (Intel): not supported
* Windows: Python 3.9-3.14


## Building from source

Please see [the Wiki compiling section](https://github.com/bertiniteam/b2/wiki/Installation) for instructions on compiling Bertini 2.

---

# Other information

The official project repository is hosted on GitHub at [github.com/bertiniteam/b2](https://github.com/bertiniteam/b2).

Please note that this is a long-term project, and is under active development.  If you want to help, please see [the wiki](https://github.com/bertiniteam/b2/wiki) for contact information.  We have opportinuties for all skill levels and interests.

---

# License

Bertini 2 is Free and Open Source Software.  Source is available under GPL Version 3, with additional terms as permitted under Section 7.


---

# Thank yous

A huge thank you to:

* HongKee Moon, for his help in getting this package to be pip-installable with a comprehensive CI build system.  
* Jack Hagen, for helping get away from the autotools and replacing CMake.
* Mike Mumm, for helping with straight-line-programs
* Jeb Collins, for writing much of the parser system, implementing the predictors, and so much more
* Tim Hodges, for contributing to the endgame implementations
* Alan Liddel, for tons of help with Python parts
* Dan Bates, Jon Hauenstein, for critical advise support, and guidance
* Andrew Sommese, Charles Wampler, Dan Bates, Jon Hauenstein, for writing Bertini 1

And to all the other people who have contributed to this package over the years.
