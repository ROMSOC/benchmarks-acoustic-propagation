# FEM Convergence Check

The scripts in this folder are to verify the consistency of the mesh, the variational form 
and the absorption of the PML.  

Verified compatibility with FEniCS 2018.0.1

## `complex.py` :
Complex container class for easier handling arithmetic between imaginary numbers and functions. 

## `planewave_sphere_classical_PML.py` :

Plane wave exact solutions are imposed on the outer domain boundary and the exterior of the object, 
while the Helmholtz equations are solved within the domain to recover the plane wave field. Errors are
computed from the exact solution to perform a convergence analysis for different mesh widths.

It is a script to verify the boundary conditions and the variational formulation with the PML.
It includes the PML domain enveloping the fluid domain, but the absorption factors are 0. 

## `monopole_sphere_classical_PML.py` :

Helmholtz equations are solved with a (-n active) finite PML region with the classical quadratic
PML absoption funcions. Monopole field is assumed around the origin placed in an non-reflecting fluid
domain ensured by the PML. Errors are computed based on the exact solution to tune the absorption strength
of the PML for different frequencies. 

To automate the optimal absorption strength, a Newton's descent algorithm is also included (although this
may not always converge to a global optimum).

## `monopole_sphere_optimal_PML.py`:

Similar script as above but with a singular PML absorption function based on Bermudez et al. (2007). 
It is optimal for all frequencies and does not need tuning. 

