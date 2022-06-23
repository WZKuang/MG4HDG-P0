# MG4HDG-P0
Numerical experiments code of multigrid preconditioners for condensd lowest-order HDG schemes with projected jumps (HDG-P0) for the reaction-diffusion equations and the generalized Stokes equations

## Requirements:
+ Netgen/NGSolve: version: 6.2.2105-289-g53df468f0, website: ngsolve.org.

## Files:
+ mymg.py: Multigrid method class.
+ diff-quad-rate.py: Multigrid preconditioned CG solver for condensed HDG-P0 for reaction-diffusion equation with know solution, convergence rate and robustness of preconditioners checked.
