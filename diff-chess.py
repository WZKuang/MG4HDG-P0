from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
# from netgen.csg import *
import matplotlib.pyplot as plt
from prol import *
from mymg import *
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver
from myIterSolver import IterSolver

# from ngsolve.webgui import Draw
# ==========coefficients==========
dim = 2
smStep = 2
gaussSm = 1
precond = 1
drawResults = True
import sys
maxLevel = int(sys.argv[1])
ratio = 1e4 # base = 1
# ===============================
smType = 'gs' if gaussSm else 'jc'
iniN = 4
if dim == 2:
    mesh = MakeStructured2DMesh(quads=False, nx=iniN, ny=iniN)
    mesh_aux = MakeStructured2DMesh(quads=False, nx=iniN * 2**maxLevel, ny=iniN * 2**maxLevel)
    maxdofs=5e7
else:
    mesh = MakeStructured3DMesh(hexes=False, nx=iniN, ny=iniN, nz=iniN)
    mesh_aux = MakeStructured2DMesh(quads=False, nx=iniN * 2**maxLevel, ny=iniN * 2**maxLevel,
                                                                         nz=iniN * 2**maxLevel)
    maxdofs=5e7
ne = mesh.ne


W0_aux = L2(mesh_aux)
gfu_aux = GridFunction(W0_aux)
for i in range(iniN * 2**maxLevel):
    for j in range(iniN * 2**maxLevel):
        loc = iniN * 2**maxLevel * i + j
        if (i+j) % 2 == 0:
            gfu_aux.vec.FV().NumPy()[2*loc:2*(loc+1)] = ratio
        else:
            gfu_aux.vec.FV().NumPy()[2*loc:2*(loc+1)] = 1
# lam_func = ratio * sin(iniN * 2**maxLevel*pi*x) * sin(iniN * 2**maxLevel*pi*y) + ratio
# gfu_aux.Set(lam_func, bonus_intorder = 4)
if drawResults:
    import netgen.gui
    Draw(gfu_aux, mesh_aux, 'var', min=0, max=ratio*1.1)
    input('continue?')
lam_avg = (gfu_aux.vec.data[1]+gfu_aux.vec.data[2]) / 2

# W0 = L2(mesh)
# gfu = GridFunction(W0)
# gfu.vec.data[:] = lam_avg
# f0 = LinearForm(W0)
# ir_c = IntegrationRule(points = [(0.1, 0.1), (0.45, 0.45), (0.2, 0.7), (0.7, 0.2)],
#                            weights= [1/8, 1/8, 1/8, 1/8])
# f0 += gfu_aux * W0.TestFunction() * dx(intrules = {TRIG:ir_c})
# f0.Assemble()
# gfu.vec.data = W0.InvM() * f0.vec
# gfu.Set(lam_func, bonus_intorder=10)
# Draw(gfu, mesh, 'L2', min=0, max=ratio*1.5)


# HDG
V = VectorL2(mesh) 
W = L2(mesh, order=1)
M = FacetFESpace(mesh, dirichlet='.*')
fes = M*V*W
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
(uhat, q, u),(vhat,r, v) = fes.TnT()  

 
# set diffusion and reaction coeff
W0 = L2(mesh)
lam_var = GridFunction(W0)
lam_var.vec.data[:] = lam_avg
# lam_var.Set(gfu_aux)
# lam_var.vec.FV().NumPy()[0::2] = 1/4 * ratio + 3/4
# lam_var.vec.FV().NumPy()[1::2] = 1/4 + 3/4 * ratio
c_rac = CoefficientFunction(1)

# true original operator
if mesh.dim == 3:
    ir_c = IntegrationRule(points = [(1/3, 1/3, 0), (1/3, 0, 1/3), (0, 1/3, 1/3), (1/3, 1/3, 1/3)],
                           weights = [1/16, 1/16, 1/16, 1/16])
    a += (1/lam_var * q * r + c_rac * u * v) * dx(intrules = {QUAD:ir_c})
elif mesh.dim == 2:
    ir_c = IntegrationRule(points = [(0.5, 0), (0, 0.5), (0.5, 0.5)],
                           weights= [1/6, 1/6, 1/6])
    a += (1/lam_var * q * r + c_rac * u * v) * dx(intrules = {TRIG:ir_c})
if mesh.dim ==2:
    ir = IntegrationRule(SEGM, 1) 
else:
    ir = IntegrationRule(TRIG, 1)
a += (-uhat*r*n+q*n*vhat
      + lam_var * alpha/h * (u-uhat)*(v-vhat)) * dx(element_boundary=True,
                                                    intrules={SEGM:ir})

# heat-source in inner subdomain
f = LinearForm(fes)
f += 1 * v * dx

a.Assemble()
# =========preconditioner init==============
pre = MultiGrid(a.mat, prol, nc=M.ndof,
                coarsedofs=fes.FreeDofs(True), w1=0.5, 
                nsmooth=smStep, sm=smType, var=False)

gfu = GridFunction(fes)
uhath, qh, uh = gfu.components

def SolveBVP(level):
    with TaskManager():
        fes.Update()
        gfu.Update()
        W0.Update()
        f.Assemble()
        if level > 0:
            et.Update()
            lam_var.Update()
            lam_var.vec.data[:] = lam_avg
            if level == maxLevel:
                lam_var.Set(gfu_aux)
#             lam_var.Set(gfu_aux)
#             lam_var.vec.FV().NumPy()[0::2] = 1
#             lam_var.vec.FV().NumPy()[1::2] = ratio
#             lam_var.vec.FV().NumPy()[2::4] = 1
#             lam_var.vec.FV().NumPy()[3::4] = 2
            
            a.Assemble()
            # point smoother
            pp = [fes.FreeDofs(True), M.ndof, [], fes.FreeDofs(True)]
            pre.Update(a.mat, pp)
        
        lams = EigenValues_Preconditioner(mat=a.mat, pre=pre)
        if precond:
            inv = CGSolver(a.mat, pre, printing=False, tol=1e-8, maxiter=100)
        else:
            inv = IterSolver(mat=a.mat, pre=pre, printrates=False, tol=1e-8, maxiter=200, 
                             freedofs=fes.FreeDofs(True))
        #inv = a.mat.Inverse(fes.FreeDofs(True))
        f.vec.data += a.harmonic_extension_trans * f.vec  
        gfu.vec.data = inv*f.vec
        it = inv.iterations   
        #it = 1
        gfu.vec.data += a.harmonic_extension * gfu.vec 
        gfu.vec.data += a.inner_solve * f.vec
        print(f"Precond:{precond}, IT: {it}, cond: {max(lams)/min(lams):.2e}, NDOFS: {M.ndof}")
        if drawResults:
        # if level == maxLevel and drawResults:
            Draw(lam_var, mesh, 'lam_var', min=0, max=ratio*1.1)
            input('continue?')
            Draw(uh, mesh, 'sol', min=0)
            input('continue?')


level = 0
print('==============================')
print(f'DIM: {dim}, Smoother Type: {smType}, Smooth Steps: {smStep}')
print('==============================')
print(f'level: {level}')
SolveBVP(level)
print('==============================')
while M.ndof < maxdofs and level < maxLevel:  
    mesh.ngmesh.Refine()
    meshRate = 2

    level += 1
    print(f'level: {level}')
    
    SolveBVP(level)
    print('==============================')