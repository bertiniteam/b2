Bertini 2:  The redevelopment of Bertini in C++.

==

## Quick links

- [Wiki](https://github.com/bertiniteam/b2/wiki)
- [Overview](#Overview)
- [Current Capabilites](#Current-capabilites)
- [Other information](#Other-information)

Thanks for checking out Bertini 2!

==

## Overview

The solution of arbitrary polynomial systems is an area of active research, and has many applications in math, science and engineering.  This program, Bertini 2, builds on the success of the first Bertini program, and seeks to eventually replace it entirely, as a powerful numerical engine.

The theoretical basis for the solution of polynomials with Bertini is "homotopy continuation", the act of "continuing" from one system into another through a "homotopy", as depicted in the below diagram.

![homotopy continuation](core/doc/images/homotopycontinuation_generic_40ppi.png "homotopy continuation")

==

## Current Capabilites

Bertini2 currently has implemented the foundations of Numerical Algebraic Geometry.  Development is ongoing, but here's what we have so far:

- C++ and Python bindings for access into any functionality.
- Construction of polynomial and non-polynomial multivariate systems.
- Evaluation of systems and Jacobians in double and arbitrary multiple precision.
- Construction of the Total Degree start system.
- Construction of homotopies (they're just systems!).
- Tracking of a start point x_0, corresponding to a particular time $t_0 \in \mathbb{C}^n$ in a homotopy $H$, from $t_0$ to $t_1$.

Development is ongoing, and we want your help!

==

## Other information

The project is hosted on [GitHub](https://github.com/bertiniteam/b2).

Please note that this is a long-term project, and is under active development.  If you want to help, please see [the wiki](https://github.com/bertiniteam/b2/wiki) for contact information.  We have opportinuties for all skill levels and interests.

Bertini 2 is Free and Open Source Software.  Source is available under GPL Version 3, with additional terms as permitted under Section 7.
