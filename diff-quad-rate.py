from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import *
import matplotlib.pyplot as plt
from prol import *
from mymg import *
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver
from myIterSolver import IterSolver

# ==========coefficients==========
import sys
if len(sys.argv) < 6:
    print('Not enough arguments: dim + sm steps + block? + GS sm? + precond?')
    exit(1)
dim = int(sys.argv[1])
smStep = int(sys.argv[2])
block = bool(int(sys.argv[3]))
gaussSm = bool(int(sys.argv[4]))
precond = bool(int(sys.argv[5]))
drawResults = False
adaptive = False
# ===============================
smType = 'gs' if gaussSm else 'jc'
if dim == 2:
    mesh = Mesh(unit_square.GenerateMesh(maxh=1/4))
    maxdofs=1e6
else:
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1/4))
    maxdofs=1e6
ne = mesh.ne


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

# set diffusion and reaction coeff to be 1
lam = CoefficientFunction(1 + 1/2 * sin (x) * sin(y))
c_rac = CoefficientFunction(1 + 1/2 * sin(x) * sin(y))

# true original operator
if mesh.dim == 3:
    ir_c = IntegrationRule(points = [(1/3, 1/3, 0), (1/3, 0, 1/3), (0, 1/3, 1/3), (1/3, 1/3, 1/3)],
                           weights = [1/16, 1/16, 1/16, 1/16])
    a += (1/lam * q * r + c_rac * u * v) * dx(intrules = {QUAD:ir_c})
elif mesh.dim == 2:
    ir_c = IntegrationRule(points = [(0.5, 0), (0, 0.5), (0.5, 0.5)],
                           weights= [1/6, 1/6, 1/6])
    a += (1/lam * q * r + c_rac * u * v) * dx(intrules = {TRIG:ir_c})
# for auxiliary operator, one-point integration
# on boundaries to replace reaction integration term over elements
if mesh.dim ==2:
    ir = IntegrationRule(SEGM, 1) 
    a += (-uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                         intrules={SEGM:ir})

else:
    ir = IntegrationRule(TRIG, 1)
    a += (-uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                          intrules={TRIG:ir})
# heat-source in inner subdomain
f = LinearForm(fes)
if dim == 2:
    f += ((1 + 1/2 * sin(x) * sin(y)) * (16 * (x - x**2) * (y - y**2)) * v
          + 32 * (y - y**2) * v
          - 8 * sin(y) * (y - y**2) * (cos(x) * (1 - 2 * x) - 2 * sin(x)) * v
          + 32 * (x - x**2) * v
          - 8 * sin(x) * (x - x**2) * (cos(y) * (1 - 2 * y) - 2 * sin(y)) * v) * dx(intrules = {TRIG:ir_c})
    # exact solution
    u_exact = 16 * x * (1 - x) * y * (1 - y)
    q_exact = CF((1 + 1/2 * sin(x) * sin(y)) * (16 * (y - y**2) * (1 - 2 * x), 16 * (x - x**2) * (1 - 2 * y)))
elif dim == 3:
    f += ((1 + 1/2 * sin(x) * sin(y)) * 16 * (x - x**2) * (y - y**2) * (z - z**2) * v
          + 32 * (y - y**2) * (z - z**2) * v
          - 8 * sin(y) * (y - y**2) * (z - z**2) * (cos(x) * (1 - 2 * x) - 2 * sin(x)) * v
          + 32 * (x - x**2) * (z - z**2) * v
          - 8 * sin(x) * (x - x**2) * (z - z**2) * (cos(y) * (1 - 2 * y) - 2 * sin(y)) * v
          + 32 * (x - x**2) * (y - y**2) * v
          + 16 * sin(x) * sin(y) * (x - x**2) * (y - y**2) * v) * dx(intrules = {QUAD:ir_c})
    u_exact = 16 * x * (1- x) * y * (1 - y) * z * (1 - z)
    q_exact = CF((1 + 1/2 * sin(x) * sin(y)) * (16 * (y - y**2) * (z - z**2) * (1 - 2 * x), 16 * (x - x**2) * (z - z**2) * (1 - 2 * y), 
                  16 * (x - x**2) * (y - y**2) * (1 - 2 * z)))

a.Assemble()
# =========preconditioner init==============
pre = MultiGrid(a.mat, prol, nc=M.ndof,
                coarsedofs=fes.FreeDofs(True), w1=0.5, 
                nsmooth=smStep, sm=smType, var=False)

gfu = GridFunction(fes)
uhath, qh, uh = gfu.components
# vtk = VTKOutput(ma=mesh,coefs=[uh], names=["sol"], filename="adpt"+str(dim),subdivision=0)

def SolveBVP():
    with TaskManager():
        fes.Update()
        gfu.Update()
        a.Assemble()
        f.Assemble()
        if mesh.ne > ne:
            et.Update()
            # ===========point or block============
            # block
            if block:
                pp = [fes.FreeDofs(True), M.ndof, [], a.mat.CreateBlockSmoother(VertexPatchBlocks(mesh,M))]
            # point
            else:
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
            import netgen.gui
            Draw(Norm(uh), mesh, 'primal norm')
            input('continue?')

l = []    # l = list of estimated total error

# postprocessing
VT = HDiv(mesh, order=0, discontinuous=True)
WT = H1(mesh)
qhs = GridFunction(VT, "flux")
uhs = GridFunction(WT, "primal")
qt, rt = VT.TnT()
at = BilinearForm(VT)
at += qt*n*rt*n*dx(element_boundary=True)
ft = LinearForm(VT)
ft += (qh*n-lam*alpha/h*(uh-uhath))*rt*n*dx(element_boundary=True, intrules={SEGM:ir})



# a posterior error estimator
def CalcError():
    with TaskManager():
        # compute the flux:
        VT.Update()      
        qhs.Update()
        at.Assemble()
        ft.Assemble()
        qhs.vec.data = at.mat.Inverse()*ft.vec
        
        WT.Update()
        uhs.Update()
        uhs.Set(uh)
        # compute estimator:
        err = 1/lam*(qh-qhs)**2 + lam*(qh/lam-grad(uhs))**2
        eta2 = Integrate(err, mesh, VOL, element_wise=True)
        maxerr = max(eta2)
        l.append ((fes.ndof, sqrt(sum(eta2))))
    #     print("ndof =", fes.ndof, " maxerr =", maxerr**0.5)
        
        # mark for refinement:
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, eta2[el.nr] > 0.25*maxerr)


level = 0
print('==============================')
print(f'DIM: {dim}, Smoother Type: {smType}, Smooth Steps: {smStep}, Block Smoother: {block}')
print('==============================')
print(f'level: {level}')
SolveBVP()
L2_uErr =  sqrt(Integrate((uh - u_exact)*(uh - u_exact), mesh))
L2_qErr =  sqrt(Integrate((qh - q_exact)*(qh - q_exact), mesh))
print(f"uh L2-error: {L2_uErr:.3E}")
print(f"qh L2-error: {L2_qErr:.3E}")
print('==============================')
prev_uErr, prev_qErr = L2_uErr, L2_qErr
u_rate, q_rate = 0, 0
while M.ndof < maxdofs:  
    if adaptive:  CalcError()
    if dim == 2:
        if adaptive:
            mesh.Refine()
        else:
            mesh.ngmesh.Refine()
        meshRate = 2
    else:
        mesh.Refine(onlyonce = True)
        meshRate = sqrt(2)

    level += 1
    SolveBVP()
    # ===== convergence check =====
    L2_uErr =  sqrt(Integrate((uh - u_exact)*(uh - u_exact), mesh))
    L2_qErr =  sqrt(Integrate((qh - q_exact)*(qh - q_exact), mesh))
    u_rate = log(prev_uErr / L2_uErr) / log(meshRate)
    q_rate = log(prev_qErr / L2_qErr) / log(meshRate)
    prev_uErr, prev_qErr = L2_uErr, L2_qErr
    print(f'level: {level}')
    print(f"uh L2-error: {L2_uErr:.3E}, uh conv rate: {u_rate:.2E}")
    print(f"qh L2-error: {L2_qErr:.3E}, qh conv rate: {q_rate:.2E}")
    print('==============================')

