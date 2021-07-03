# Planewave scattering by rigid sphere

An incident plane wave is impinged on a rigid sphere, and these scripts computes the scattered field
utilizing a classical quadratic PML absorption boundary condition. The variational form is written 
in terms of the scattered field. Exact solution is known and the errors can be computed to validate
the solver.

The absorption coeffients for the quadratic PML functions are estimated using the monopole test case 
(`source/01_fem_convergence_check/monopole_sphere_classical_PML`) and used here.

Exact Solution is provided by C++ scripts provided in the `exact` folder which are compiled during 
runtime by the FEniCS provided JIT compiler. 

FE Errors can be computed and convergence behavior can be investigated using these scripts.