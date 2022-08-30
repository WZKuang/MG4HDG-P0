from ngsolve.meshes import MakeStructured3DMesh
from ngsolve import *
from netgen.csg import unit_cube
import numpy as np
from ngsolve.krylovspace import CGSolver
import time as timeit
from ngsolve.la import EigenValues_Preconditioner
from numpy import pi
###### import prolongation
ngsglobals.msg_level=0

#SetNumThreads (8)

def HDivTest(structure=False, h0 = 0.3,
           tol = 1e-10, atol = 1e-14, nlvls =6,
           c_vis=1, c_div=1, c_lo=1e-4,
           printing=False, alpha=4, draw=False, adapt = False):
    print("#########################################################################")
    print("c_vis: ", c_vis, "c_div: ", c_div, "c_lo: ", c_lo)
    print("adaptive: ", adapt)
    print("#########################################################################")

    with TaskManager():
         if structure:
             mesh = MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2) 
         else:
             mesh = Mesh(unit_cube.GenerateMesh(maxh=h0))
         # reference solution
         uex0 = -2*sin(pi*x)*cos(pi*y)*cos(pi*z)
         uex1 = cos(pi*x)*sin(pi*y)*cos(pi*z)
         uex2 = cos(pi*x)*cos(pi*y)*sin(pi*z)
         uex = CoefficientFunction((uex0, uex1, uex2))
         lapuex0 = (uex0.Diff(x)).Diff(x)+ (uex0.Diff(y)).Diff(y)\
                 +(uex0.Diff(z)).Diff(z)
         lapuex1 = (uex1.Diff(x)).Diff(x)+ (uex1.Diff(y)).Diff(y) \
                 +(uex1.Diff(z)).Diff(z)
         lapuex2 = (uex2.Diff(x)).Diff(x)+ (uex2.Diff(y)).Diff(y) \
                 +(uex2.Diff(z)).Diff(z)
         source = CoefficientFunction((-c_vis/2*lapuex0+c_lo*uex0, 
             -c_vis/2*lapuex1+c_lo*uex1, 
             -c_vis/2*lapuex2+c_lo*uex2))

         V = FESpace("BDM1", mesh, dirichlet=".*")
         M = FESpace("HCurlP1", mesh)
         fes = FESpace([V,M])
             
         gfu = GridFunction(fes)
         uh, uhath = gfu.components

         space_flux = HDiv(mesh, order= 1)
         gf_flux0 = GridFunction(space_flux)
         gf_flux1 = GridFunction(space_flux)
         gf_flux2 = GridFunction(space_flux)

         maxits = 2000
         # the bilinear-form
         (u,uhat), (v,vhat) = fes.TnT()
         a = BilinearForm(fes, symmetric=True)
         
         # viscous part
         def tang(v):
             return v-(v*n)*n
         h = specialcf.mesh_size
         n = specialcf.normal(mesh.dim)
         # viscous term (laplacian)
         a += c_vis*InnerProduct(Sym(Grad(u)), Sym(Grad(v)))*dx
         a += (c_div*div(u)*div(v)+c_lo*u*v)*dx # Hdiv-elliptic problem
         jmpu = tang(u-uhat)
         jmpv = tang(v-vhat)
         a += c_vis*(-Sym(Grad(u))*n*jmpv-Sym(Grad(v))*n*jmpu
                 +alpha/h*jmpu*jmpv)*dx(element_boundary=True)
         
         c = Preconditioner(a, "multigrid", smoother="block", 
                 #coarsetype="direct",
                 #coarsesmoothingsteps=1, 
                 #updatealways=False,
                 #cycle=2,
                 smoothingsteps=3,
                 #test=True, 
                 inverse="umfpack"
                 #inverse="sparsecholesky"
                 )
         # linear form
         f = LinearForm(fes)
         
         f += source*v* dx

         a.Assemble()
         
         def CalcZZError():
             # compute the flux:
             space_flux.Update()      
             gf_flux0.Update()
             gf_flux1.Update()
             gf_flux2.Update()
             dgfu = Grad(gfu)
             flux0 =  CoefficientFunction((dgfu[0], dgfu[1], dgfu[2]))
             flux1 =  CoefficientFunction((dgfu[3], dgfu[4], dgfu[5]))
             flux2 =  CoefficientFunction((dgfu[6], dgfu[7], dgfu[8]))
             gf_flux0.Set(flux0)
             gf_flux1.Set(flux1)
             gf_flux2.Set(flux2)
             # compute estimator:
             err = ((flux0-gf_flux0)**2+(flux1-gf_flux1)**2
                     +(flux2-gf_flux2)**2)

             eta2 = Integrate(err, mesh, VOL, element_wise=True)
             maxerr = max(eta2)
             # mark for refinement (TYPE 2)
             sumE = sum(eta2)
             markerr = 0
             nmark = 0
             # set-false
             for el in mesh.Elements():    
               mesh.SetRefinementFlag(el, False)
             for el in mesh.Elements(BND):    
               mesh.SetRefinementFlag(el, False)
             while markerr < 0.25*sumE:
               for el in mesh.Elements():
                 if eta2[el.nr] > 0.81*maxerr :
                   mesh.SetRefinementFlag(el, True)
                   markerr += eta2[el.nr]
                   nmark += 1
                   eta2[el.nr] = 0.
               maxerr = max(eta2)
             return nmark # return # of marked elements


         # solve soln on coarse mesh
         a.Assemble()
         f.Assemble()
         gfu.vec.data = a.mat.Inverse(fes.FreeDofs())*f.vec
         
         if draw:
             Draw(uh, mesh, "uh")
             Draw(uhath, mesh, "uhath")
             input("??")
         
         errU0, errUhat0 = 1, 1
         count = 0
         
         # coarse grid op.
         #for k in range(nlvls):
         while fes.ndof < 1e6: # coarse grid operators
             t0 = timeit.time()
             if adapt:
                 nmark = CalcZZError()
             else:
                 nmark = mesh.ne
             mesh.Refine()
             fes.Update()
             gfu.Update()
             a.Assemble()
             f.Assemble()
             t1 = timeit.time()
             # eigenvalues
             lams = EigenValues_Preconditioner(mat=a.mat, pre=c.mat)
             inv = CGSolver(a.mat, c.mat, printing=printing, tol=tol, 
                 abstol=atol, maxsteps=maxits)
             t2 = timeit.time()
             gfu.vec.data = inv*f.vec
             # condensation part
             t3 = timeit.time()
             errU = sqrt(Integrate((uh-uex)**2*dx, mesh))
             errUhat = sqrt(Integrate(tang(uhath-uex)**2/h
                  *dx(element_boundary=True), mesh))
             rateU = -log(errU/errU0)/log(2)
             rateUhat = -log(errUhat/errUhat0)/log(2)
             errU0, errUhat0 = errU, errUhat
             print("asmb: %.2e eigen: %.2e solve: %.2e  iter: %5i ndof: %.2e nele: %.2e"%(
                         t1-t0, t2-t1, t3-t2,inv.iterations, fes.ndof,
                         mesh.ne))
             print("errU: %.2e rateU: %.2f errUhat: %.2e  rateUhat:%.2f lamMax: %.3e lamMin: %.3e cond: %.3e"%(
                         errU, rateU, errUhat, rateUhat, 
                         max(lams), min(lams),max(lams)/(min(lams)+1e-14)))

             if inv.iterations==maxits:
               print("SOLVER FAILED")
               break
             if draw:
                 gf = GridFunction(VectorH1(mesh, order=2))
                 gf.Set(gfu)
                 Draw(gf, mesh, "X")
                 input("??")


HDivTest(structure=False, h0 = 0.3, alpha=6, c_vis=1, c_div=1e4, c_lo=0, printing=True)
