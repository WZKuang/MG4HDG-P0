# MG4HDG-P0
Numerical experiments code of multigrid (MG) preconditioners for condensd lowest-order HDG schemes with projected jumps (HDG-P0) for the reaction-diffusion equations and the generalized Stokes equations

Link to the paper: https://arxiv.org/abs/2208.14418
## Requirements:
+ **Netgen/NGSolve**: version: 6.2.2105-289-g53df468f0, website: ngsolve.org.

## Files:
+ **prol**: NGSolve add on for intergrid transfer of hybrid facet space, needs to be built and installed.
+ **mymg.py**: Multigrid method class.
+ **myIterSolver.py**: Iteration solver class.
+ **diff-quad-rate.py**: Solve HDG-P0 for reaction-diffusion equations with known solution.
+ **diff-jump.py**: Solve HDG-P0 for reaction-diffusion equations with jump coefficients on the coarsest mesh.
+ **diff-chess.py**: Solve HDG-P0 for reaction-diffusion equations with jump coefficients on the finest mesh in 2D.
+ **stokes-quad-rate.py**: Solve HDG-P0 for the generalized Stokes equations with known solution.
+ **stokes-lid.py**: One-time augmented Lagrangian Uzawa iteration to solve the lid-driven cavity problems.
+ **stokes-back.py**: One-time augmented Lagrangian Uzawa iteration to solve the backward-facing step flow problem.