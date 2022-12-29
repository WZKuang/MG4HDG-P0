from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import *
from prol import *
from mymg import *
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver

import time as timeit

# ============ parameters ===========
c_vis, c_div= 1, 1e8

import sys
dim = int(sys.argv[1]) #2
ns = int(sys.argv[2]) #1 #number of smoothersa
c_lo = int(sys.argv[3])
var= int(sys.argv[4]) #True
wc = int(sys.argv[5]) #False
# ==================================
adapt = True
drawResults = True
js = False  # THIS IS A HACK # Nope #

# === geometric setup === #
L = 5
maxLevel = 25
if dim == 2:
    maxdofs = 3e6
    c0 = 0.3826
    ## Backwards-facing step flow geo
    geo = SplineGeometry()
    pnts = [(0.5, 0), (L, 0), (L, 1), (0, 1), (0, 0.5), (0.5, 0.5)]
    pind = [geo.AppendPoint(*pnt) for pnt in pnts]
    geo.Append(['line', pind[0], pind[1]], leftdomain=1, rightdomain=0, bc="wall")
    geo.Append(['line', pind[1], pind[2]], leftdomain=1, rightdomain=0, bc="outlet")
    geo.Append(['line', pind[2], pind[3]], leftdomain=1, rightdomain=0, bc="wall")
    geo.Append(['line', pind[3], pind[4]], leftdomain=1, rightdomain=0, bc="inlet")
    geo.Append(['line', pind[4], pind[5]], leftdomain=1, rightdomain=0, bc="wall")
    geo.Append(['line', pind[5], pind[0]], leftdomain=1, rightdomain=0, bc="wall")
    
    mesh = Mesh(geo.GenerateMesh(maxh=1/2))
    uin = CoefficientFunction(16 * (1 - y) * (y - 0.5))
else:
    maxdofs = 1e6
    c0 = 0.241
    ## Backwards-facing step flow geo
    geo = CSGeometry()

    left = Plane(Pnt(0, 0, 0), Vec(-1, 0, 0)).bc("inlet")
    right = Plane(Pnt(L, 0, 0), Vec(1, 0, 0)).bc("outlet")
    box1 = OrthoBrick(Pnt(-1, 0, 0), Pnt(L + 1, 1, 1)).bc("wall")
    box2 = OrthoBrick(Pnt(-1, -1, -1), Pnt(0.5, 3, 0.5)).bc("wall")

    geo.Add(box1 * left * right - box2)
    mesh = Mesh(geo.GenerateMesh(maxh=1/1.8))
    uin = CoefficientFunction(64 * (1 - z) * (z - 0.5) * (1 - y) * y)
ne = mesh.ne

# HDG
V = MatrixValued(L2(mesh), mesh.dim, False)
W = VectorL2(mesh, order=1)
M = FacetFESpace(mesh, dirichlet='wall|inlet')
Q = L2(mesh)
if mesh.dim == 2:
    fes = M * M * V * W * Q
else:
    fes = M * M * M * V * W * Q

et = meshTopology(mesh, mesh.dim)
et.Update()
if mesh.dim == 2:
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
    ir = IntegrationRule(TRIG, 1)
    a += (-uhat * (G * n) + (L * n) * vhat
          + uhat * n * q - p * vhat * n
          + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                              intrules={TRIG: ir})
    a_ax += (c_lo * h / 3 * uhat * vhat
             - uhat * (G * n) + (L * n) * vhat
             + uhat * n * q - p * vhat * n
             + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                                 intrules={TRIG: ir})

f = LinearForm(fes)
a.Assemble()
# a_ax.Assemble()
# jacobi smoother
pre = MultiGrid(a.mat, prol, nc=M.ndof,
                coarsedofs=fes.FreeDofs(True), w1=0.8,
                nsmooth=ns, sm="gs", var=var,
                he=True, dim=dim, wcycle=wc, js=js)

gfu = GridFunction(fes)
if mesh.dim == 2:
    uhath0, uhath1, Lh, uh, ph = gfu.components
    uhath = CF((uhath0, uhath1))
else:
    uhath0, uhath1, uhath2, Lh, uh, ph = gfu.components
    uhath = CF((uhath0, uhath1, uhath2))


# vtk = VTKOutput(ma=mesh,coefs=[uh, ph], names=["vel", "pres"], 
#         filename="Stokes"+str(dim)+str(adapt),subdivision=0)

def SolveBVP(level):
    with TaskManager():
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
                pdofs[j * M.ndof:(j + 1) * M.ndof] = inner
            # he_prol
            pp.append(a.mat.Inverse(pdofs, inverse="sparsecholesky"))
            # bk smoother
            bjac = et.CreateSmoother(a, {"blocktype": "vertexpatch"})
            pp.append(bjac)
            pre.Update(a.mat, pp)
        t2 = timeit.time()
        # estimate condition number
        lams = EigenValues_Preconditioner(mat=a.mat, pre=pre)
        #lams = np.array([1,1])
        t21 = timeit.time()
        inv = CGSolver(a.mat, pre, printing=False, tol=1e-8, maxiter=40)
        # === dirichlet BC === #
        if dim == 2:
            gfu.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))
        else:
            gfu.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))

        f.vec.data -= a.mat * gfu.vec
        f.vec.data += a.harmonic_extension_trans * f.vec
        
        gfu.vec.data += inv * f.vec
        it = inv.iterations
        #gfu.vec.data += a.mat.Inverse(fes.FreeDofs(True))*f.vec
        #it = 1
        
        gfu.vec.data += a.harmonic_extension * gfu.vec
        gfu.vec.data += a.inner_solve * f.vec
        t3 = timeit.time()

        # print(f"Time to find EIG Val: {t21 - t2}, MAX PREC LAM: {max(lams)}; MIN PREC LAM: {min(lams)}")
        print("IT: %2.0i" % (it), "cond: %.2e" % (max(lams) / min(lams)), "NDOFS: %.2e" % (
            sum(fes.FreeDofs(True))),
            "Assemble: %.2e Prec: %.2e Solve: %.2e" % (t1 - t0, t2 - t1, t3 - t21))
        
        if drawResults:
            import netgen.gui
            Draw(Norm(uh), mesh, 'velocity norm')
            input('continue?')


l = []  # l = list of estimated total error

# # postprocessing
VT0 = HDiv(mesh, order=0, discontinuous=True)
if dim == 2:
    VT = VT0 * VT0
else:
    VT = VT0 * VT0 * VT0

WT = VectorH1(mesh)
Lh0 = GridFunction(VT, "flux")
uhs = GridFunction(WT, "primal")
qts, rts = VT.TnT()
if dim == 2:
    qt = CF((qts[0], qts[1]), dims=(2, 2))
    rt = CF((rts[0], rts[1]), dims=(2, 2))
    Lhs = CF((Lh0.components[0], Lh0.components[1]), dims=(2, 2))
else:
    qt = CF((qts[0], qts[1], qts[2]), dims=(3, 3))
    rt = CF((rts[0], rts[1], rts[2]), dims=(3, 3))
    Lhs = CF((Lh0.components[0], Lh0.components[1], Lh0.components[2]), dims=(3, 3))

at = BilinearForm(VT)
ft = LinearForm(VT)
at += (qt * n) * (rt * n) * dx(element_boundary=True)
ft += (Lh * n - c_vis * alpha / h * (uh - uhath)) * (rt * n) * dx(element_boundary=True, intrules={SEGM: ir})


# # a posterior error estimator
def CalcError():
    with TaskManager():
        # compute the flux:
        VT.Update()
        Lh0.Update()
        at.Assemble()
        ft.Assemble()
        Lh0.vec.data = at.mat.Inverse(inverse="sparsecholesky") * ft.vec

        WT.Update()
        uhs.Update()
        uhs.Set(uh)
        # compute estimator:
        err = 1 / c_vis * InnerProduct(Lh - Lhs, Lh - Lhs) + c_vis * InnerProduct(Lh / c_vis - Grad(uhs),
                                                                                Lh / c_vis - Grad(uhs))
        # err += 1 / c0 ** 2 * div(uhs) ** 2
        eta2 = Integrate(err, mesh, VOL, element_wise=True)
        maxerr = max(eta2)
        l.append((fes.ndof, sqrt(sum(eta2))))
        #     print("ndof =", fes.ndof, " maxerr =", maxerr**0.5)

        # mark for refinement:
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, eta2[el.nr] > 0.75 * maxerr)


print('=====================')
print(f"DIM: {dim}, Adapt: {adapt}, vertex-patch-GS steps: {ns}, var-V: {var}, W-cycle: {wc}, c_low: {c_lo:.1E}, eps: {c_div:.1E}")
print('=====================') 
SolveBVP(0)
level = 0
while level < maxLevel and M.ndof * dim < maxdofs:
    level += 1
    print('=========')
    print(f'Level: {level}')
    if adapt == True:
        CalcError()
    if dim == 3:
        mesh.Refine(onlyonce=True)
        # if adapt:
        #     mesh.Refine()
        # else:
        #     mesh.Refine(onlyonce=True)
    elif dim == 2:
        if adapt:  mesh.Refine()
        else: mesh.ngmesh.Refine()
    # exit if total global dofs exceed a tol 
    SolveBVP(level)
