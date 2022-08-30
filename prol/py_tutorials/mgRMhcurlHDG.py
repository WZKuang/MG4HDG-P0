######## NOTE: DG does not work nicely on unstructured mesh (???)
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from netgen.csg import * 
from netgen.geom2d import unit_square
import numpy as np
from ngsolve.krylovspace import CGSolver
import time as timeit
###### import prolongation
from prol import * 
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
             mesh = unit_square.GenerateMesh(maxh=h0)
             mesh.EnableTable("parentedges")
             mesh = Mesh(mesh)
             #mesh = MakeStructured2DMesh(quads=False, nx= nx, ny=nx)
             #V = HCurl(mesh, order= order, dirichlet=".*")
             V = FESpace("HCurlP1X", mesh, dirichlet=".*")
             W = H1(mesh, order = order+1, dirichlet=".*")
             #M = HDiv(mesh, order = order-1, dirichlet=".*") # hybrid XX
             M = HDiv(mesh, order = order, dirichlet=".*") # hybrid XX
             fes = FESpace([V,W, M])
             x1, y1 = x*(x-1), y*(y-1)
             x2, y2 = 2*x-1, 2*y-1
             xx, yy = 5*x**2-5*x+1, 5*y**2-5*y+1
             scale = 4**6
             uex0 = scale*y1**3*x1**2*x2
             uex1 = scale*x1**3*y1**2*y2
             uex = CoefficientFunction((uex0, uex1))
             wex = scale*(1/3*x1**3*y1**3 - 2*t**2/5/(1-nu)*(
                     y1**3*x1*xx+x1**3*y1*yy))
             dwex = CoefficientFunction((wex.Diff(x), wex.Diff(y)))
             source = scale*E/(1-nu**2)*(
                     y1*xx*(2*y1**2+x1*yy)
                    +x1*yy*(2*x1**2+y1*xx))
         gfu = GridFunction(fes)
         uh, wh, uhath = gfu.components
             

         mt = meshTopology(mesh, mesh.dim)
         mt.Update() # always update
         #prol1 = EdgeProlongationTrig(mesh, mt, True)
         prol1 = V.Prolongation() 
         prol2 = NodalProlongationTrig(mesh, mt, True)
         prol3 = EdgeProlongationTrig2(mesh, mt, True)
         # one-point integration
         ir = IntegrationRule(SEGM, 2*order+1) # reduced integration

         maxits = 4000
         # the bilinear-form
         (theta, w, that), (eta, v, ehat) = fes.TnT()
         a = BilinearForm(fes, symmetric=True)
         
         h = specialcf.mesh_size
         n = specialcf.normal(mesh.dim)
         
         # Normal velocity
         def norm(v):
             return (v*n)*n
         
         # stress in terms of eps
         def sigma(eps):
             return fac*((1-nu)*eps+nu*Trace(eps)*Id(mesh.dim))

         # viscous term (laplacian)
         eps = Sym(Grad(theta))
         stress = sigma(eps)
         
         deps = Sym(Grad(eta))
         dstress = sigma(deps)
         a += InnerProduct(stress, Grad(eta))*dx

         jmpu = norm(theta-that)
         jmpv = norm(eta-ehat)

         a += (-stress*n*jmpv-dstress*n*jmpu
                 +alpha/h*jmpu*jmpv)*dx(element_boundary=True,intrules={SEGM:ir})
         a += lam*(grad(w)-theta)*(grad(v)-eta)*dx # shear energy
         
         # linear form
         f = LinearForm(fes)
         
         if test ==1:
             f += source*v * dx

         a.Assemble()
         pre = MultiGrid3(a.mat, prol1, prol2, prol3, uh.vec.size, wh.vec.size,
                 uhath.vec.size, coarsedofs=fes.FreeDofs(),
             nsmooth=ns, wcycle=wcycle, var=var, sm=sm, w1=weight1)
         

         # solve soln on coarse mesh
         a.Assemble()
         f.Assemble()
         gfu.vec.data = a.mat.Inverse(fes.FreeDofs())*f.vec
         
         if draw:
             Draw(wex, mesh, "x0")
             Draw(uex, mesh, "y0")
             Draw(wh, mesh, "x")
             Draw(uh, mesh, "y")
             input("??")
         
         errW0, errW20, errU0 = 1, 1, 1
         count = 0
         # coarse grid op.
         #for k in range(nlvls):
         while fes.ndof < 4e6: # coarse grid operators
             t0 = timeit.time()
             count += 1
             nmark = mesh.ne
             mesh.Refine()
             fes.Update()
             gfu.Update()
             mt.Update()
             a.Assemble()
             t1 = timeit.time()
             
             # blocks (THIS IS EXPENSIVE) 
             if bk=="vp": # vertex patch (costly)
                 eblocks = VertexPatchBlocks(mesh, fes)
             elif bk=="ep":
                 eblocks = EdgePatchBlocks(mesh, fes)
             pre.Update(a.mat, eblocks, uh.vec.size, wh.vec.size, 
                     uhath.vec.size)
             
             ##### XXX
             t2 = timeit.time()
             
             t3 = timeit.time()
             f.Assemble()
             inv = CGSolver(a.mat, pre, printing=printing, tol=tol, 
                 abstol=atol, maxsteps=maxits)
             gfu.vec.data = inv*f.vec
             # condensation part
             t4 = timeit.time()
             errW = sqrt(Integrate((wh-wex)**2, mesh))
             errW2 = sqrt(Integrate((grad(wh)-dwex)**2, mesh))
             errU = sqrt(Integrate((uh-uex)**2, mesh))
             rateW = -log(errW/errW0)/log(2)
             rateW2 = -log(errW2/errW20)/log(2)
             rateU = -log(errU/errU0)/log(2)
             #print("errU:%.4e rate: %.3f errW2:%.4e rate: %.3f errW:%.4e rate: %.3f "%(
             #    sqrt(errU),rateU, sqrt(errW2), rateW2,  sqrt(errW), rateW))
             errW0, errW20, errU0 = errW, errW2, errU
             print("assemble: %.2e blocking: %.2e solve: %.2e iter: %5i ndof: %.2e nele: %.2e"%(
                         t1-t0, t2-t1, t4-t3,inv.iterations, fes.ndof,
                         mesh.ne))
             if inv.iterations==maxits:
               print("SOLVER FAILED")
               break
             if draw:
                 Redraw()
                 input("??")


HCurlTest(
        test=1, order=1, h0 = 1, nx=2, alpha=4, structure=False,slip=True,
        bk = "vp", wcycle=0, var=False, sm="gs", ns=4, weight1=0.1,
        t=1e-0, E=12, nu=0, printing=False, draw=False, 
        adapt=False)
