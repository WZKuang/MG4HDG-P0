from ngsolve.meshes import MakeStructured3DMesh
from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import numpy as np
from ngsolve.krylovspace import CGSolver
import time as timeit
from ngsolve.la import EigenValues_Preconditioner
from numpy import pi
###### import prolongation
ngsglobals.msg_level=4

def HDivTest(test= 1,structure=False, nx=10, h0 = 0.3,
           tol = 1e-10, atol = 1e-14, nlvls =6,
           c_vis=1, c_div=1, c_lo=1e-4,
           printing=False, alpha=4, draw=False, adapt = False,
           blocktype="edgepatch", cycle=1, sm=1, ism = 1):
    print("#########################################################################")
    print("c_vis: ", c_vis, "c_div: ", c_div, "c_lo: ", c_lo)
    print("blocktype: ", blocktype, " cycle: ", cycle, " sm: ", sm, " ism: ", ism)
    print("structured: ", structure)
    print("#########################################################################")

    with TaskManager(pajetrace=10**8):
         if test == 11 : #2D 
             mesh = Mesh(unit_square.GenerateMesh(maxh=h0))
         else:
             if structure:
                 mesh = MakeStructured3DMesh(hexes=False, nx=nx, ny=nx, nz=nx) 
             else:
                 mesh = Mesh(unit_cube.GenerateMesh(maxh=h0))
         SetTestoutFile ("test"+str(test)+"DG"+blocktype+".out")
         if test==1:
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
         elif test ==2: # cavity
            # dirichlet data 
            uex = CoefficientFunction((1, 0, 0))
         elif test==11:
            # reference solution
            uex0 = -sin(pi*x)*cos(pi*y)
            uex1 = cos(pi*x)*sin(pi*y)
            uex = CoefficientFunction((uex0, uex1))
            lapuex0 = (uex0.Diff(x)).Diff(x)+ (uex0.Diff(y)).Diff(y)
            lapuex1 = (uex1.Diff(x)).Diff(x)+ (uex1.Diff(y)).Diff(y)
            source = CoefficientFunction((-c_vis/2*lapuex0+c_lo*uex0, 
                -c_vis/2*lapuex1+c_lo*uex1)) 

         if c_vis>0:
           fes = FESpace("BDM1", mesh, dirichlet=".*", dgjumps=True)
         else:
           fes = FESpace("BDM1", mesh, dirichlet=".*")

         gfu = GridFunction(fes)

         space_flux = HDiv(mesh, order= 1)
         gf_flux0 = GridFunction(space_flux)
         gf_flux1 = GridFunction(space_flux)
         gf_flux2 = GridFunction(space_flux)

         maxits = 1000
         # the bilinear-form
         u, v = fes.TnT()
         a = BilinearForm(fes, symmetric=True)
         a += (c_div*div(u)*div(v)+c_lo*u*v)*dx # Hdiv-elliptic problem
         # viscous part
         def tang(v):
             return v-(v*n)*n
         h = specialcf.mesh_size
         n = specialcf.normal(mesh.dim)
         # viscous term (laplacian)
         if c_vis>0:
            ir = IntegrationRule(TRIG, 1)
            a += c_vis*InnerProduct(Sym(Grad(u)), Sym(Grad(v)))*dx
            avgdu = 0.5*(Sym(Grad(u))+Sym(Grad(u.Other())))*n
            avgdv = 0.5*(Sym(Grad(v))+Sym(Grad(v.Other())))*n
            jmpu = tang(u-u.Other())
            jmpv = tang(v-v.Other())
            stab = alpha*(1/h+1/h.Other())
            a += c_vis*(-avgdu*jmpv-avgdv*jmpu
                    +stab*jmpu*jmpv)*dx(skeleton=True)
         # linear form
         f = LinearForm(fes)
         if test==1 or test==11:
           f += source*v* dx
         elif test==2:
           stab0 = 2*alpha/h
           a += c_vis*(-Sym(Grad(u))*n*tang(v)-Sym(Grad(v))*n*tang(u)
                    +stab0*tang(u)*tang(v))*ds(skeleton=True)
           f += c_vis*(-Sym(Grad(v))*n*tang(uex)
                    +stab0*tang(uex)*tang(v))*ds(skeleton=True,
                        definedon=mesh.Boundaries("top"))
         c = Preconditioner(a, "multigrid", smoother="block", 
                 blocktype=blocktype,
                 cycle=cycle,
                 smoothingsteps=sm,
                 increasesmoothingsteps=ism,
                 test=True, 
                 #inverse="umfpack"
                 inverse="sparsecholesky"
                 )
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
             Draw(gfu, mesh, "uh")
             input("??")
         
         errU0, errUhat0 = 1, 1
         count = 0
         ne0 = mesh.ne
         # coarse grid op.
         #for k in range(nlvls):
         while fes.ndof < 1e6: # coarse grid operators
             t0 = timeit.time()
             if adapt:
                 nmark = CalcZZError()
             else:
                 nmark = mesh.ne
             #mesh.Refine(onlyonce=True)
             mesh.Refine()
             fes.Update()
             gfu.Update()
             a.Assemble()
             #continue # HAHAHA
             f.Assemble()
             t1 = timeit.time()
             inv = CGSolver(a.mat, c.mat, printing=printing, tol=tol, 
                 abstol=atol, maxsteps=maxits)
             t2 = timeit.time()
             gfu.vec.data = inv*f.vec
             # condensation part
             t3 = timeit.time()
             errU = sqrt(Integrate((gfu-uex)**2*dx, mesh))
             rateU = -log(errU/errU0)/log(mesh.ne/ne0)*mesh.dim
             ne0 = mesh.ne
             errU0 = errU
             print("asmb: %.2e eigen: %.2e solve: %.2e  iter: %5i ndof: %.2e nele: %.2e"%(
                         t1-t0, t2-t1, t3-t2,inv.iterations, fes.ndof,
                         mesh.ne), 
                         "errU: %.2e rateU: %.2f"%(
                         errU, rateU))

             if inv.iterations==maxits:
               print("SOLVER FAILED")
               break
             if draw:
               Redraw()  
               input("??")

#structure = [True, False]
#coefs = [(0,1,1), (0, 1e4, 1), (1, 1, 0), (1, 1e4, 0)]
#bks = ["edgepatch", "vertexpatch"]
structure = [True]
#coefs = [(1, 1, 0), (1, 1e4, 0)]
coefs = [(1, 1e4, 0)]
#coefs = [(0, 1e4, 1)]
bks = ["vertexpatch"]
#bks = ["edgepatch", "vertexpatch"]
tests = [11]
for it in tests:
  for st in structure:
    for coef in coefs:
      for bk in bks:
        HDivTest(test=it, structure=st, nx = 10, h0 = 0.5, alpha=6, 
          c_vis=coef[0], c_div=coef[1], c_lo=coef[2],
          blocktype=bk, cycle=1, sm=1, ism=1, adapt=False, printing=False,
          draw=False)
