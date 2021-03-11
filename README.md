# 2d-cft-single-strip-joining-quench
Code for calculations about a local quench in 2d CFT involving taking two ends of a strip and joining them.

This code accompanies section 4.3 in [`arXiv:2103.`](https://arxiv.org/2103.), but should be helpful in a more general context.

Calculations in this quench set-up involves, on the first hand, calculating a specific conformal transformation between an annulus and a second 'quench' manifold, known as a doubly connected Schwarz-Christoffel map.[1]
The quench manifold is an plane with two slits on the real line, arranged to be symmetric under an inversion about a circle of radius `1`; the slits nearly touch each other.
Since the width of the annulus is a conformal invariant, it is a property of the quench manifold that has to be calculated.

This typically has to be done numerically, and we use a previously written package DSCPACK. [2]

A short description of the various files:
1. `ang1em2driver.f` is a fortran code that uses DSCPACK to invert the map on the time-reflection-symmetric slice for a quench manifold in which the ends of the two slits are a distance `.01` from each other. It can be easily modified to do the same for different distances.
2. `schwdriver.f` calculates the derivatives of the conformal map for a quench manifold in which the distance between the slits is `1E-5`. It saves the derivative of the map between the annulus and the quench manifold above to `wprods.txt`; and also the derivatives of a closely related map from a finite cylinder to an infinite cylinder with reflection-symmetric slits to `sprime.txt`
3. `calcschw.py` uses formulae from [1] to calculate the Schwarzian of the conformal map.

----

[1]: DeLillo, T. K., Elcrat, A. R., & Pfaltzgraff, J. A. (2001). Schwarz--Christoffel Mapping of the Annulus. SIAM review, 43(3), 469-477. <https://doi.org/10.1137/S0036144500375280>.  
[2]: Chenglie Hu. 1998. Algorithm 785: a software package for computing Schwarz-Christoffel conformal transformation for doubly connected polygonal regions. ACM Trans. Math. Softw. 24, 3 (Sept. 1998), 317–333. <https://doi.org/10.1145/292395.291204>.  
Code: <http://www.netlib.org/toms-2014-06-10/785>.
