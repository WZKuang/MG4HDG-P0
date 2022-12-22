from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import *
import matplotlib.pyplot as plt
from prol import *
from mymg import *
#from mymgHack import *
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver, MinResSolver, GMResSolver

import time as timeit
# ============ parameters ===========
c_vis, c_div= 1, 1e8

import sys
dim = int(sys.argv[1]) #2
ns = int(sys.argv[2]) #1 #number of smoothersa
c_lo = int(sys.argv[3])
var= int(sys.argv[4]) #True
wc = int(sys.argv[5]) #False

adapt = False
js = False   # THIS IS A HACK # Nope #
# ==================================

if dim == 2:
    mesh = Mesh(unit_square.GenerateMesh(maxh=1/4))
    maxdofs=1e6
    c0 = 0.3826
else:
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1/2))
    maxdofs=1e6
    c0 = 0.241
ne = mesh.ne

# HDG
V = MatrixValued(L2(mesh),mesh.dim, False)
W = VectorL2(mesh, order=1)
M = FacetFESpace(mesh, dirichlet=".*")
Q = L2(mesh)
if mesh.dim==2: 
    fes = M*M*V*W*Q
else:
    fes = M*M*M*V*W*Q

et = meshTopology(mesh, mesh.dim)
et.Update() 
if mesh.dim==2:
    prol = FacetProlongationTrig(mesh, et) 
else:
    prol = FacetProlongationTet(mesh, et) 

n = specialcf.normal(mesh.dim) 
h = specialcf.mesh_size
alpha = 1
a = BilinearForm(fes, symmetric=False, condense=True)
a_ax = BilinearForm(fes, symmetric=False, condense=True)

if mesh.dim == 2:
    (uhat0, uhat1, L, u, p), (vhat0, vhat1, G, v, q) = fes.TnT()
    uhat = CF((uhat0, uhat1))
    vhat = CF((vhat0, vhat1))
    ir_c = IntegrationRule(points = [(0.5, 0), (0, 0.5), (0.5, 0.5)],
                           weights= [1/6, 1/6, 1/6])
    a += (1/c_vis * InnerProduct(L,G) + c_lo * u * v + 1/c_div * p * q) * dx(intrules = {TRIG:ir_c})
    a_ax += (1/c_vis * InnerProduct(L,G) + 1/c_div * p * q) * dx(intrules = {TRIG:ir_c})
    ir = IntegrationRule(SEGM, 1)
    a += (-uhat * (G * n) + (L * n) * vhat
          + uhat * n * q - p * vhat * n
          + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                              intrules={SEGM: ir})

    a_ax += (c_lo * h / 3 * uhat * vhat
             - uhat * (G * n) + (L * n) * vhat
             + uhat * n * q - p * vhat * n
             + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                                 intrules={SEGM: ir})
elif mesh.dim == 3:
    (uhat0, uhat1, uhat2, L, u, p), (vhat0, vhat1, vhat2, G, v, q) = fes.TnT()
    uhat = CF((uhat0, uhat1, uhat2))
    vhat = CF((vhat0, vhat1, vhat2))
    ir = IntegrationRule(TRIG, 1)
    ir_c = IntegrationRule(points=[(1 / 3, 1 / 3, 0), (1 / 3, 0, 1 / 3), (0, 1 / 3, 1 / 3), (1 / 3, 1 / 3, 1 / 3)],
                           weights=[1 / 16, 1 / 16, 1 / 16, 1 / 16])
    a += (1 / c_vis * InnerProduct(L, G) + c_lo * u * v + 1 / c_div * p * q) * dx(intrules={QUAD: ir_c})
    a_ax += (1 / c_vis * InnerProduct(L, G) + 1 / c_div * p * q) * dx(intrules={QUAD: ir_c})
    a += (-uhat*(G*n)+(L*n)*vhat
          + uhat*n*q - p*vhat*n
          + c_vis*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                  intrules={TRIG:ir})

    a_ax += (c_lo * h / 3 * uhat * vhat
             - uhat*(G*n)+(L*n)*vhat
             + uhat*n*q - p*vhat*n
             + c_vis*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                  intrules={TRIG:ir})

f = LinearForm(fes)
a.Assemble()
# a_ax.Assemble()
# jacobi smoother
pre = MultiGrid(a.mat, prol, nc=M.ndof,
                coarsedofs=fes.FreeDofs(True), w1=0.8, 
                nsmooth=ns, sm="gs", var=var, 
               he=True, dim=dim, wcycle=wc, js=js)

gfu = GridFunction(fes)
if mesh.dim==2:
    uhath0, uhath1, Lh, uh, ph = gfu.components
    uhath = CF((uhath0, uhath1))
else:
    uhath0, uhath1, uhath2, Lh, uh, ph = gfu.components
    uhath = CF((uhath0, uhath1, uhath2))
    
# vtk = VTKOutput(ma=mesh,coefs=[uh, ph], names=["vel", "pres"], 
#         filename="Stokes"+str(dim)+str(adapt),subdivision=0)

def SolveBVP(level):
    t0 = timeit.time()
    fes.Update()
    gfu.Update()
    a.Assemble()
    # a_ax.Assemble()
    f.Assemble()
    t1 = timeit.time()
    
    # print(f"Global DOFs: {sum(fes.FreeDofs(True))}, Vhat Dofs: {sum(M.FreeDofs()) * dim}")
    ## Update MG
    if level > 0:
        et.Update()
        pp = [fes.FreeDofs(True), M.ndof]
        pdofs = BitArray(fes.ndof) 
        pdofs[:] = 0 
        inner = prol.GetInnerDofs(level)
        for j in range(dim):
            pdofs[j*M.ndof:(j+1)*M.ndof] = inner
        # he_prol
        pp.append(a.mat.Inverse(pdofs, inverse="sparsecholesky"))
        # bk smoother
        bjac = et.CreateSmoother(a, {"blocktype":"vertexpatch"})
        pp.append(bjac)
        pre.Update(a.mat, pp)
    t2 = timeit.time()
    # estimate condition number
    lams = EigenValues_Preconditioner(mat=a.mat, pre=pre)
    #lams = np.array([1,1])
    t21 = timeit.time()
    inv = CGSolver(a.mat, pre, printing=False, tol=1e-8, maxiter=60)
    # dirichlet BC
    if dim ==2:
        gfu.components[0].Set(4*x*(1-x), definedon=mesh.Boundaries("top"))
    else:
        gfu.components[0].Set(16*x*(1-x)*y*(1-y), definedon=mesh.Boundaries("top"))
    
    f.vec.data -= a.mat*gfu.vec    
    f.vec.data += a.harmonic_extension_trans * f.vec  
    gfu.vec.data += inv*f.vec
    it = inv.iterations
#     gfu.vec.data += a.mat.Inverse(fes.FreeDofs(True))*f.vec
#     it = 1
    gfu.vec.data += a.harmonic_extension * gfu.vec 
    gfu.vec.data += a.inner_solve * f.vec
    t3 = timeit.time()

    # print(f"Time to find EIG Val: {t21 - t2}, MAX PREC LAM: {max(lams)}; MIN PREC LAM: {min(lams)}")
    print("IT: %2.0i"%(it), "cond: %.2e"%(max(lams)/min(lams)), "NDOFS: %.2e"%(
        sum(fes.FreeDofs(True))), 
          "Assemble: %.2e Prec: %.2e Solve: %.2e"%(t1-t0, t2-t1, t3-t21))

# SolveBVP(0)
level = 1
while True:
    with TaskManager():
        if adapt==True:
            mesh.Refine()
        else:
            if dim == 3:
                mesh.Refine(onlyonce = True)
            else:
                mesh.ngmesh.Refine()
        # exit if total global dofs exceed a tol 
        M.Update()
        if (M.ndof*dim > maxdofs):
            #print(M.ndof*dim)
            break
        if level==1:
            print('=====================')
            print(f"DIM: {dim}, Adapt: {adapt}, vertex-patch-GS steps: {ns}, var-V: {var}, W-cycle: {wc}, c_low: {c_lo:.1E}, eps: {c_div:.1E}")
            print('=====================')
        SolveBVP(level)
        # from ngsolve.webgui import Draw
        # import netgen.gui
        # Draw(uh, mesh, 'norm')
        # input('continue?')
        level +=1




