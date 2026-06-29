
# Quick links

- [Documentation](https://bertini2.org)
- [`bertini2` PyPI package](https://pypi.org/project/bertini2/)
- [Wiki](https://github.com/bertiniteam/b2/wiki)

---

# Overview

The solution of arbitrary polynomial systems is an area of active research, and has many applications in math, science and engineering.  This software, Bertini 2, is a re-implementation of [Bertini 1](https://bertini.nd.edu) from C into C++/Python.

The theoretical basis for the solution of polynomials with Bertini is a theorem which bounds the number of solutions a system may have. It sits together with the numerical computational tool of "homotopy continuation". the act of "continuing" from one system into another through a "homotopy", as depicted in the below diagram:

<img src="https://raw.githubusercontent.com/bertiniteam/b2/develop/doc_resources/images/homotopycontinuation_generic.png" alt="homotopy continuation"
	title="homotopy continuation" width="400" height="285"/>

---

# Current capabilites

Bertini 2 has already implemented most of the foundations of Numerical Algebraic Geometry.  Development is ongoing, and here's what we have so far:

- A blackbox that implements tracktype 0 (zerodim solving) for the total degree and mhom start systems.  User homotopy not yet exposed
- Through Python bindings, runtime scriptable construction of systems and interactivity with their zero-dimensional solutions.
- Python construction of multivariate polynomial and non-polynomial systems, including from linear algebra operations via numpy and Bertini 2's variables and extended numeric types.  
- Evaluation of systems and their Jacobians in double and arbitrary multiple precision.
- Construction of the Total Degree and Multihomogeneous start systems.
- Construction of homotopies (they're just systems with path variables defined).
- Tracking of a start point x_0, corresponding to a particular time $t_0 \in \mathbb{C}^n$ in a homotopy $H$, from $t_0$ to $t_1$.
- Running of the Power Series and Cauchy endgames, in double, multiple, and adaptive precision.

See the published docs for latest the officially released documentation (which will always lag behind the doc from the dev version).

---

# Missing functionality

✨ In-progress!!!  

* Numerical irreducible decomposition
* Membership testing
* and other algorithms

Users wanting a more stable implementation are recommended to use [Bertini 1](https://bertini.nd.edu) or [homotopycontinuation.jl](https://www.juliahomotopycontinuation.org/), or one of the other packages implementing the theory.  

---

# Installation

## Pre-built wheels -- the way to go!

The Python package `bertini2` provides pre-built wheels for Linux, macOS, and Windows. 


```bash
pip install bertini2
```

Once it's installed, you `import bertini` (not `import bertini2`!)

Installing the wheel also puts the classic blackbox command-line solver on your PATH as
`bertini2` (threads-only -- run `bertini2 --help`).  For MPI / multi-rank cluster
parallelism, build the CLI from source.

* Linux: Python 3.10-3.14
* MacOS (Apple Silicon): Python 3.10-3.14
* MacOS (Intel): not currently supported
* Windows: Python 3.10-3.14


## Building from source

Please see [the Wiki compiling section](https://github.com/bertiniteam/b2/wiki/Installation) for instructions on compiling Bertini 2.

## The `bertini2` command-line solver

The classic blackbox CLI (`bertini2 input`, with a `CONFIG`/`INPUT` file) is a **pure C++**
program -- it does not depend on Python -- and every build path delivers it:

| You want... | Do this | You get |
| --- | --- | --- |
| The easy way (threads-only) | `pip install bertini2` | `bertini2` on PATH + the `bertini` Python package |
| From source, with Python | `pip install .` (in a clone) | same as the wheel, built locally |
| From source, with Python, editable | `pip install -e .` | `bertini2` resolves to the in-tree build |
| **Just the solver, no Python** | `cmake -B build -S . && cmake --build build --target bertini2_exe && cmake --install build --component cli --prefix <prefix>` | only `<prefix>/bin/bertini2` |
| Full C++ dev install | `cmake --install build` | `bin/bertini2` + `libbertini2` + headers |

A plain `cmake` build skips the Python bindings (`BUILD_PYTHON_BINDINGS` is off by default), so the
CLI-only path never needs a Python interpreter. The binary statically links our code but still needs
the usual shared libraries on the system (Boost, GMP, MPFR, MPC, and MPI for multi-rank runs) -- the
wheel vendors these; a from-source build expects them present (e.g. via your HPC modules or a conda
env). The wheel CLI is **threads-only**; for MPI / multi-rank cluster parallelism, build from source.

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
