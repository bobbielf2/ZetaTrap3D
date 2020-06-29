# High-order locally corrected trapezoidal rule for Laplace and Helmholtz layer potentials on surfaces in 3D

This is the MATLAB code accompanying the paper: 

* B. Wu and P.G. Martinsson, Corrected Trapezoidal Rules for Boundary Integral Equations in Three Dimensions (2020, to be published)

Let's consider a Laplace or Helmholtz layer potential from a smooth surface to a target location on the same surface, then the kernel is singular at the target location. Assume that the surface is parameterized over a rectangular domain and discretized uniformly using the double trapezoidal rule. Our corrected trapezoidal rule approximates this on-surface potential to high accuracy, provided that either the surface parameterization is doubly-periodic (such as a toroidal surface) or the potential is compactly supported on the surface.

Author: Bowei Wu, 2020/6/28

Also contain supporting functions modified from Alex Barnett's [BIE3D](https://github.com/ahbarnett/BIE3D) package

### Description of the main test files:

* `test_laplace3d_on_quartic_patch.m` : evaluation of the Laplace layer potentials compactly supported on a surface patch. This routine contains a self convergence test.
* `test_helmholtz3d_on_quartic_patch.m` evaluation of the Helmholtz layer potentials compactly supported on a surface patch. This routine contains a self convergence test.
* `test_laplace3d_bie.m` : solution of the Laplace Dirichlet and Neumann problems exterior to a toroidal surface. This routine contains a convergence test.
* `test_helmholtz3d_bie.m`  : solution of the Helmholtz Dirichlet and Neumann problems exterior to a toroidal surface. This routine contains a convergence test.

Supporting functions:

* `epstein_zeta.m` : evaluation of the Epstein zeta function and its parametric derivatives.
* `Lap3dLocCorr.m` : construct matrices associated with the Laplace layer potentials on a surface using locally corrected trapezoidal rule
* `Helm3dLocCorr.m` : construct matrices associated with the Helmholtz layer potentials on a surface using locally corrected trapezoidal rule
* `Lap3dSLPmat.m`, `Lap3dDLPmat.m` : construct matrices associated with Laplace layer potentials using native quadratures
* `Helm3dSLPmat.m`, `Helm3dDLPmat.m` : construct matrices associated with Helmholtz layer potentials using native quadratures
* `quadr_doubleptr.m`, `wobblytorus.m`, `showsurf.m` are (modified) functions from [BIE3D](https://github.com/ahbarnett/BIE3D)

