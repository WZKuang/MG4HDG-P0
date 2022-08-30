######## NOTE: HCurl 2D
from ngsolve import *
from meshes import MakeStructured2DMesh, MakeStructured3DMesh
from netgen.csg import * 
from netgen.geom2d import unit_square
import numpy as np
from ngsolve.krylovspace import CGSolver
import time as timeit
###### import prolongation
from mymg import * # multigrid stuff
import netgen.gui
ngsglobals.msg_level=0
SetHeapSize(int(4e9))

def HCurlTest(test=1, nx=2, order=1, structure=False, 
           ns=1, bk="vp", tol = 1e-10, atol = 1e-14,
           wcycle=0, var=False, sm="gs", weight1=0.1, nlvls =6,
           t=1e-2, E=1, nu=0.4, h0 = 0.3, slip=True,
           printing=False, alpha=4, alphaA=False, draw=False, adapt = False):
    print("#########################################################################")
    print("test: ", test, "order: ", order)
    print("nu: ", nu, "t: ", t)
    print("wcycle: ", wcycle, " var: ", var, " sm: ", sm, " w1: ", weight1)
    print("Blocks: ", bk, "; #smooth: ", ns, "; CG tol(rel err): ", tol)
    print("adaptive: ", adapt)
    print("#########################################################################")

    with TaskManager():
         # parameters
         lam = 5*E/12/(1+nu)*t**(-2)
         fac = E/12/(1-nu**2)
         alpha *= fac # scale stabilization parameter

         if test ==1: # clamped plate
             mesh = unit_square.GenerateMesh(maxh=1)
             mesh.EnableTable ("parentedges")
             mesh = Mesh(mesh)
             fes = FESpace("HCurlP1", mesh, dirichlet=".*", dgjumps=True)
         
         gfu = GridFunction(fes)
         
         # P1-Hcurl prolongation
         prol1 = fes.Prolongation()

         # one-point integration
         ir = IntegrationRule(SEGM, 2*order-1) # reduced integration

         maxits = 100
         # the bilinear-form
         theta, eta = fes.TnT()
         a = BilinearForm(fes, symmetric=True)
         
         h = specialcf.mesh_size
         n = specialcf.normal(mesh.dim)
         
         # Normal velocity
         def norm(v):
             return (v*n)*n

         # viscous term (laplacian)
         a += InnerProduct(Sym(Grad(theta)), Sym(Grad(eta)))*dx
         avgdu = 0.5*(Sym(Grad(theta))+Sym(Grad(theta.Other())))*n
         avgdv = 0.5*(Sym(Grad(eta))+Sym(Grad(eta.Other())))*n
         jmpu = norm(theta-theta.Other())
         jmpv = norm(eta-eta.Other())
         a += (-avgdu*jmpv-avgdv*jmpu+alpha/h*jmpu*jmpv)*dx(skeleton=True)
         a += (1e4*curl(theta)*curl(eta))*dx
         
         # linear form
         f = LinearForm(fes)
         
         f += eta[0] * dx

         a.Assemble()
         pre = MultiGrid1(a.mat, prol1, fes.ndof,
                 coarsedofs=fes.FreeDofs(),
             nsmooth=ns, wcycle=wcycle, var=var, sm=sm, w1=weight1)
         

         # solve soln on coarse mesh
         a.Assemble()
         f.Assemble()
         gfu.vec.data = a.mat.Inverse(fes.FreeDofs())*f.vec
         
         errW0, errW20, errU0 = 1, 1, 1
         count = 0
         # coarse grid op.
         #for k in range(nlvls):
         while fes.ndof < 1e5: # coarse grid operators
             t0 = timeit.time()
             count += 1
             nmark = mesh.ne
             mesh.Refine()
             fes.Update()
             gfu.Update()
             a.Assemble()
             t1 = timeit.time()
             
             # blocks (THIS IS EXPENSIVE) 
             if bk=="vp": # vertex patch (costly)
                 eblocks = VertexPatchBlocks(mesh, fes)
             elif bk=="ep":
                 eblocks = EdgePatchBlocks(mesh, fes)
             pre.Update(a.mat, eblocks, fes.ndof)
             
             ##### XXX
             t2 = timeit.time()
             
             t3 = timeit.time()
             f.Assemble()
             inv = CGSolver(a.mat, pre, printing=printing, tol=tol, 
                 maxiter=maxits)
             gfu.vec.data = inv*f.vec
             # condensation part
             t4 = timeit.time()
             print("assemble: %.2e blocking: %.2e solve: %.2e iter: %5i ndof: %.2e nele: %.2e"%(
                         t1-t0, t2-t1, t4-t3,inv.iterations, fes.ndof,
                         mesh.ne))
             if inv.iterations==maxits:
               print("SOLVER FAILED")
               break


HCurlTest(
        test=1, order=1, h0 = 1, nx=2, alpha=4, structure=False,slip=True,
        bk = "vp", wcycle=0, var=False, sm="gs", ns=1, weight1=0.1,
        t=1, E=12, nu=0, printing=False, draw=False, 
        adapt=False)
