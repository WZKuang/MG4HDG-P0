# MG4HDG-P0
Numerical experiments code of multigrid (MG) preconditioners for condensd lowest-order HDG schemes with projected jumps (HDG-P0) for the reaction-diffusion equations and the generalized Stokes equations

## Requirements:
+ **Netgen/NGSolve**: version: 6.2.2105-289-g53df468f0, website: ngsolve.org.

## Files:
+ **mymg.py**: Multigrid method class, head file.
+ **diff-quad-rate.py**: MG preconditioned CG solver for condensed HDG-P0 for reaction-diffusion equations with known solution.
+ **diff-Chip.py**: MG preconditioned CG solver for condensed HDG-P0 for reaction-diffusion equations with with jump coefficients.
+ **stokes-quad-rate.py**: MG preconditioned CG solver for the augmented Lagrangian Uzawa iteration method to solve the condensed HDG-P0 for generalized Stokes equations with known solution.
+ **stokes-lid.py**: MG preconditioned CG solver, one-time augmented Lagrangian Uzawa iteration to solve the lid-driven cavity problems.
+ **stokes-back.py**: MG preconditioned CG solver, one-time augmented Lagrangian Uzawa iteration to solve the backward-facing step flow problems.