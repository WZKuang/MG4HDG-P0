from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import *
from prol import *
from mymg import *
from myIterSolver import IterSolver
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver, MinResSolver, GMResSolver
import sys
import time as timeit

# ============ parameters ===========
c_vis, c_div = 1, 1e4

if len(sys.argv) < 7:
    print('Not enough arguments: dim + sm steps + c_low + var-V? + W-cycle? + precond?')
    exit(1)
dim = int(sys.argv[1])
ns = int(sys.argv[2]) #1 #number of smoothers
c_lo = float(sys.argv[3])
var = bool(int(sys.argv[4]))
wc = bool(int(sys.argv[5]))
precond = bool(int(sys.argv[6]))

uzawa = True  # whether to use AL-Uzawa or saddle point prob solver for aux operator
uzIt = 3
directSol = False
# ==================================

if dim == 2:
    mesh = Mesh(unit_square.GenerateMesh(maxh=1 / 4))
    maxdofs = 1e5
elif dim == 3:
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1 / 4))
    maxdofs = 2e7

ne = mesh.ne

# HDG
V = MatrixValued(L2(mesh), mesh.dim, False)
W = VectorL2(mesh, order=1)
M = FacetFESpace(mesh, dirichlet=".*")
Q = L2(mesh, lowest_order_wb=True)  # global pressure in original prob, no stabilization
Q0 = L2(mesh, lowest_order_wb=not uzawa)  # auxiliary operator, with stabilization

if dim == 2:
    fes = M * M * Q * V * W
    fes0 = M * M * Q0 * V * W
elif dim == 3:
    fes = M * M * M * Q * V * W
    fes0 = M * M * M * Q0 * V * W

et = meshTopology(mesh, mesh.dim)
et.Update()
if mesh.dim == 2:
    prolM = FacetProlongationTrig(mesh, et)  # Facet prol
else:
    prolM = FacetProlongationTet(mesh, et)
prolP = ElementProlongation0(mesh, et)  # constant pressure prol

n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size
alpha = 1
# beta = 1 / (1 + c_lo * h**2 / c_vis / (dim + 1))
a = BilinearForm(fes, symmetric=False, condense=True)
f = LinearForm(fes)
# auxiliary operator, provide SPD AL-Uzawa solver or saddle MG solver
a0 = BilinearForm(fes0, symmetric=False, condense=True)
f0 = LinearForm(fes0)
if mesh.dim == 2:
    ir_c = IntegrationRule(points=[(0.5, 0), (0, 0.5), (0.5, 0.5)],
                           weights=[1 / 6, 1 / 6, 1 / 6])
    (uhatx, uhaty, p, L, u), (vhatx, vhaty, q, G, v) = fes.TnT()
    uhat = CF((uhatx, uhaty))
    vhat = CF((vhatx, vhaty))
    a += (1 / c_vis * InnerProduct(L, G) + c_lo * u * v) * dx(intrules={TRIG: ir_c})
    ir = IntegrationRule(SEGM, 1)
    a += (-uhat * (G * n) + (L * n) * vhat
          - uhat * n * q - p * vhat * n
          + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                              intrules={SEGM: ir})

    (uhatx0, uhaty0, p0, L0, u0), (vhatx0, vhaty0, q0, G0, v0) = fes0.TnT()
    uhat0 = CF((uhatx0, uhaty0))
    vhat0 = CF((vhatx0, vhaty0))
    a0 += (1 / c_vis * InnerProduct(L0, G0) + c_lo * u0 * v0 - 1 / c_div * p0 * q0) * dx(intrules={TRIG: ir_c})
    ir = IntegrationRule(SEGM, 1)
    a0 += (-uhat0 * (G0 * n) + (L0 * n) * vhat0
           - uhat0 * n * q0 - p0 * vhat0 * n
           + c_vis * alpha / h * (u0 - uhat0) * (v0 - vhat0)) * dx(element_boundary=True,
                                                                   intrules={SEGM: ir})
elif mesh.dim == 3:
    ir_c = IntegrationRule(points=[(1 / 3, 1 / 3, 0), (1 / 3, 0, 1 / 3), (0, 1 / 3, 1 / 3), (1 / 3, 1 / 3, 1 / 3)],
                           weights=[1 / 16, 1 / 16, 1 / 16, 1 / 16])
    (uhatx, uhaty, uhatz, p, L, u), (vhatx, vhaty, vhatz, q, G, v) = fes.TnT()
    uhat = CF((uhatx, uhaty, uhatz))
    vhat = CF((vhatx, vhaty, vhatz))
    a += (1 / c_vis * InnerProduct(L, G) + c_lo * u * v) * dx(intrules={QUAD: ir_c})
    ir = IntegrationRule(TRIG, 1)
    a += (-uhat * (G * n) + (L * n) * vhat
          - uhat * n * q - p * vhat * n
          + c_vis * alpha / h * (u - uhat) * (v - vhat)) * dx(element_boundary=True,
                                                              intrules={TRIG: ir})

    (uhatx0, uhaty0, uhatz0, p0, L0, u0), (vhatx0, vhaty0, vhatz0, q0, G0, v0) = fes0.TnT()
    uhat0 = CF((uhatx0, uhaty0, uhatz0))
    vhat0 = CF((vhatx0, vhaty0, vhatz0))
    a0 += (1 / c_vis * InnerProduct(L0, G0) + c_lo * u0 * v0 - 1 / c_div * p0 * q0) * dx(intrules={QUAD: ir_c})
    ir = IntegrationRule(TRIG, 1)
    a0 += (-uhat0 * (G0 * n) + (L0 * n) * vhat0
           - uhat0 * n * q0 - p0 * vhat0 * n
           + c_vis * alpha / h * (u0 - uhat0) * (v0 - vhat0)) * dx(element_boundary=True,
                                                                   intrules={TRIG: ir})

gfu = GridFunction(fes)
gfu0 = GridFunction(fes0)
if mesh.dim == 2:
    uhathx, uhathy, ph, Lh, uh= gfu.components
    uhath = CF((uhathx, uhathy))

    uhathx0, uhathy0, ph0, Lh0, uh0 = gfu0.components
    uhath0 = CF((uhathx0, uhathy0))
elif mesh.dim == 3:
    uhathx, uhathy, uhathz, ph, Lh, uh = gfu.components
    uhath = CF((uhathx, uhathy, uhathz))

    uhathx0, uhathy0, uhathz0, ph0, Lh0, uh0 = gfu0.components
    uhath0 = CF((uhathx0, uhathy0, uhathz0))

# ===== exact solution
if mesh.dim == 2:
    u_exact1 = x ** 2 * (x - 1) ** 2 * 2 * y * (1 - y) * (2 * y - 1)
    u_exact2 = y ** 2 * (y - 1) ** 2 * 2 * x * (x - 1) * (2 * x - 1)
    u_exact = CF((u_exact1, u_exact2))

    L_exactXX = c_vis * (2 * x * (x - 1) ** 2 * 2 * y * (1 - y) * (2 * y - 1)
                         + x ** 2 * 2 * (x - 1) * 2 * y * (1 - y) * (2 * y - 1))
    L_exactXY = c_vis * (x ** 2 * (x - 1) ** 2 * 2 * (1 - y) * (2 * y - 1)
                         - x ** 2 * (x - 1) ** 2 * 2 * y * (2 * y - 1)
                         + 2 * x ** 2 * (x - 1) ** 2 * 2 * y * (1 - y))
    L_exactYX = c_vis * (y ** 2 * (y - 1) ** 2 * 2 * (x - 1) * (2 * x - 1)
                         + y ** 2 * (y - 1) ** 2 * 2 * x * (2 * x - 1)
                         + 2 * y ** 2 * (y - 1) ** 2 * 2 * x * (x - 1))
    L_exactYY = c_vis * (2 * y * (y - 1) ** 2 * 2 * x * (x - 1) * (2 * x - 1)
                         + y ** 2 * 2 * (y - 1) * 2 * x * (x - 1) * (2 * x - 1))
    L_exact = CF((L_exactXX, L_exactXY, L_exactYX, L_exactYY), dims=(2, 2))

    p_exact = x * (1 - x) * (1 - y) - 1 / 12

    f += (-c_vis * (4 * y * (1 - y) * (2 * y - 1) * ((1 - 2 * x) ** 2 - 2 * x * (1 - x))
                    + 12 * x ** 2 * (1 - x) ** 2 * (1 - 2 * y))
          + (1 - 2 * x) * (1 - y)) * v[0] * dx(intrules={TRIG: ir_c})
    f += (-c_vis * (4 * x * (1 - x) * (1 - 2 * x) * ((1 - 2 * y) ** 2 - 2 * y * (1 - y))
                    + 12 * y ** 2 * (1 - y) ** 2 * (2 * x - 1))
          - x * (1 - x)) * v[1] * dx(intrules={TRIG: ir_c})
    f += c_lo * u_exact * v * dx(intrules={TRIG: ir_c})

    f0 += (-c_vis * (4 * y * (1 - y) * (2 * y - 1) * ((1 - 2 * x) ** 2 - 2 * x * (1 - x))
                    + 12 * x ** 2 * (1 - x) ** 2 * (1 - 2 * y))
          + (1 - 2 * x) * (1 - y)) * v0[0] * dx(intrules={TRIG: ir_c})
    f0 += (-c_vis * (4 * x * (1 - x) * (1 - 2 * x) * ((1 - 2 * y) ** 2 - 2 * y * (1 - y))
                    + 12 * y ** 2 * (1 - y) ** 2 * (2 * x - 1))
          - x * (1 - x)) * v0[1] * dx(intrules={TRIG: ir_c})
    f0 += c_lo * u_exact * v0 * dx(intrules={TRIG: ir_c})
elif mesh.dim == 3:
    u_exact1 = x ** 2 * (x - 1) ** 2 * (2 * y - 6 * y ** 2 + 4 * y ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
    u_exact2 = y ** 2 * (y - 1) ** 2 * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
    u_exact3 = -2 * z ** 2 * (z - 1) ** 2 * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
    u_exact = CF((u_exact1, u_exact2, u_exact3))

    L_exactXX = c_vis * (2 * x * (x - 1) ** 2 * (2 * y - 6 * y ** 2 + 4 * y ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                         + x ** 2 * 2 * (x - 1) * (2 * y - 6 * y ** 2 + 4 * y ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3))
    L_exactXY = c_vis * (x ** 2 * (x - 1) ** 2 * (2 - 12 * y + 12 * y ** 2) * (2 * z - 6 * z ** 2 + 4 * z ** 3))
    L_exactXZ = c_vis * (x ** 2 * (x - 1) ** 2 * (2 - 12 * z + 12 * z ** 2) * (2 * y - 6 * y ** 2 + 4 * y ** 3))
    L_exactYX = c_vis * (y ** 2 * (y - 1) ** 2 * (2 - 12 * x + 12 * x ** 2) * (2 * z - 6 * z ** 2 + 4 * z ** 3))
    L_exactYY = c_vis * (2 * y * (y - 1) ** 2 * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                         + y ** 2 * 2 * (y - 1) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3))
    L_exactYZ = c_vis * (y ** 2 * (y - 1) ** 2 * (2 - 12 * z + 12 * z ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
    L_exactZX = c_vis * (-2 * z ** 2 * (z - 1) ** 2 * (2 - 12 * x + 12 * x ** 2) * (2 * y - 6 * y ** 2 + 4 * y ** 3))
    L_exactZY = c_vis * (-2 * z ** 2 * (z - 1) ** 2 * (2 - 12 * y + 12 * y ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
    L_exactZZ = c_vis * (
                -2 * 2 * z * (z - 1) ** 2 * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
                - 2 * z ** 2 * 2 * (z - 1) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * y - 6 * y ** 2 + 4 * y ** 3))
    L_exact = CF((L_exactXX, L_exactXY, L_exactXZ,
                  L_exactYX, L_exactYY, L_exactYZ,
                  L_exactZX, L_exactZY, L_exactZZ), dims=(3, 3))

    p_exact = x * (1 - x) * (1 - y) * (1 - z) - 1 / 24
    f += (-c_vis * ((2 - 12 * x + 12 * x ** 2) * (2 * y - 6 * y ** 2 + 4 * y ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (x ** 2 - 2 * x ** 3 + x ** 4) * (-12 + 24 * y) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (x ** 2 - 2 * x ** 3 + x ** 4) * (-12 + 24 * z) * (2 * y - 6 * y ** 2 + 4 * y ** 3))
          + (1 - 2 * x) * (1 - y) * (1 - z)
          ) * v[0] * dx(intrules={QUAD: ir_c})
    f += (-c_vis * ((2 - 12 * y + 12 * y ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (y ** 2 - 2 * y ** 3 + y ** 4) * (-12 + 24 * x) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (y ** 2 - 2 * y ** 3 + y ** 4) * (-12 + 24 * z) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
          - x * (1 - x) * (1 - z)
          ) * v[1] * dx(intrules={QUAD: ir_c})
    f += (2 * c_vis * (
                (2 - 12 * z + 12 * z ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
                + (z ** 2 - 2 * z ** 3 + z ** 4) * (-12 + 24 * x) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
                + (z ** 2 - 2 * z ** 3 + z ** 4) * (-12 + 24 * y) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
          - x * (1 - x) * (1 - y)
          ) * v[2] * dx(intrules={QUAD: ir_c})
    f += c_lo * u_exact * v * dx(intrules={QUAD: ir_c})

    f0 += (-c_vis * ((2 - 12 * x + 12 * x ** 2) * (2 * y - 6 * y ** 2 + 4 * y ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (x ** 2 - 2 * x ** 3 + x ** 4) * (-12 + 24 * y) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (x ** 2 - 2 * x ** 3 + x ** 4) * (-12 + 24 * z) * (2 * y - 6 * y ** 2 + 4 * y ** 3))
          + (1 - 2 * x) * (1 - y) * (1 - z)
          ) * v0[0] * dx(intrules={QUAD: ir_c})
    f0 += (-c_vis * ((2 - 12 * y + 12 * y ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (y ** 2 - 2 * y ** 3 + y ** 4) * (-12 + 24 * x) * (2 * z - 6 * z ** 2 + 4 * z ** 3)
                    + (y ** 2 - 2 * y ** 3 + y ** 4) * (-12 + 24 * z) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
          - x * (1 - x) * (1 - z)
          ) * v0[1] * dx(intrules={QUAD: ir_c})
    f0 += (2 * c_vis * (
            (2 - 12 * z + 12 * z ** 2) * (2 * x - 6 * x ** 2 + 4 * x ** 3) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
            + (z ** 2 - 2 * z ** 3 + z ** 4) * (-12 + 24 * x) * (2 * y - 6 * y ** 2 + 4 * y ** 3)
            + (z ** 2 - 2 * z ** 3 + z ** 4) * (-12 + 24 * y) * (2 * x - 6 * x ** 2 + 4 * x ** 3))
          - x * (1 - x) * (1 - y)
          ) * v0[2] * dx(intrules={QUAD: ir_c})
    f0 += c_lo * u_exact * v0 * dx(intrules={QUAD: ir_c})

# ==== mixed bilinear form for div and grad rec matrices
# Note: NC FESpace has different DOFs ordering
if dim == 2:
    MTmp = M * M
    (uhatTmpx, uhatTmpy), (vhatTmpx, vhatTmpy) = MTmp.TnT()
    uhatTmp = CF((uhatTmpx, uhatTmpy))
    vhatTmp = CF((vhatTmpx, vhatTmpy))
elif dim == 3:
    MTmp = M * M * M
    (uhatTmpx, uhatTmpy, uhatTmpz), (vhatTmpx, vhatTmpy, vhatTmpz) = MTmp.TnT()
    uhatTmp = CF((uhatTmpx, uhatTmpy, uhatTmpz))
    vhatTmp = CF((vhatTmpx, vhatTmpy, vhatTmpz))

pNC, qNC = Q.TnT()
# grad mat: p -> M * dim
p_M_mix = BilinearForm(trialspace = Q, testspace = MTmp)
p_M_mix += -pNC * n * vhatTmp * dx(element_boundary=True)
# p mass mat
pMass = BilinearForm(Q)
pMass += pNC * qNC * dx
# embedding mat: M * dim -> fes
MFesPdt = BilinearForm(trialspace=MTmp, testspace=fes)
MFesPdt += uhatTmp * vhat * dx(element_boundary=True)
fesMass = BilinearForm(fes)
fesMass += (u * v + p * q + InnerProduct(L, G)) * dx(intrules={TRIG: ir_c})
fesMass += (uhat * vhat) * dx(element_boundary=True)


a0.Assemble()
pFesEmb = Embedding(fes.ndof, fes.Range(2)) if dim == 2 else \
          Embedding(fes.ndof, fes.Range(3))
pMass.Assemble(); pMass_inv = pMass.mat.CreateSmoother()
p_M_mix.Assemble()
fesMass.Assemble(); MFesPdt.Assemble()
MFesEmb = fesMass.mat.CreateSmoother() @ MFesPdt.mat
pTransf = c_div * pFesEmb @ pMass_inv @ p_M_mix.mat.T @ MFesEmb.T
Mdof0 = BitArray(fes0.ndof)
Mdof0[:] = 0
Mdof0[: dim * M.ndof] = 1
he0 = [a0.mat.Inverse(Mdof0, inverse="umfpack"), pTransf]
# initialize MG preconditioner
if uzawa:
    pre = MultiGrid(a0.mat, prolM, nc=M.ndof,
                    coarsedofs=fes0.FreeDofs(True), w1=0.8,
                    nsmooth=ns, sm="gs", var=var,
                    he=True, dim=dim, wcycle=wc, he0=he0)
else:
    pre = MultiGrid(a0.mat, [prolM, prolP], nc=[M.ndof, Q0.ndof],
                    coarsedofs=fes0.FreeDofs(True), w1=0.8,
                    nsmooth=ns, sm="gs", var=var,
                    he=True, dim=dim, wcycle=wc, he0=he0)

def SolveBVP(level):
    with TaskManager():
        t0 = timeit.time()
        fes.Update(); fes0.Update()
        gfu.Update(); gfu0.Update()
        a.Assemble(); a0.Assemble()
        f.Assemble(); f0.Assemble()

        Q.Update(); MTmp.Update()
        p_M_mix.Assemble()
        pMass.Assemble()
        pMass_inv = pMass.mat.CreateSmoother()
        MFesPdt.Assemble(); fesMass.Assemble()
        fesMass_inv = fesMass.mat.CreateSmoother()
        # embedding mat: M * dim -> fes
        MFesEmb = fesMass_inv @ MFesPdt.mat
        # embedding mat: p -> fes
        pFesEmb = Embedding(fes.ndof, fes.Range(2)) if dim == 2 else Embedding(fes.ndof, fes.Range(3))

        t1 = timeit.time()

        ## Update MG
        # TODO: fix MG update for the condensed mixed formulation
        # TODO: update Q_l^T in mixed MG-he, the bug in current setting is wrong uHe
        if level > 0:
            et.Update()
            if uzawa:
                pp = [fes0.FreeDofs(True)]
                pp.append(M.ndof)
                pdofs = BitArray(fes0.ndof)
                pdofs[:] = 0
                inner = prolM.GetInnerDofs(level)
                for j in range(dim):
                    pdofs[j * M.ndof:(j + 1) * M.ndof] = inner
                # he_prol
                pp.append(a0.mat.Inverse(pdofs, inverse="sparsecholesky"))
                # bk smoother
                # bjac = et.CreateSmoother(a0, {"blocktype": "vertexpatch"})
                pp.append(fes0.CreateSmoothBlocks(vertex=True, globalDofs=True))
                pre.Update(a0.mat, pp)
            else:
                pp = [fes0.FreeDofs(True)]
                pp.append([M.ndof, Q0.ndof])
                pdofs = BitArray(fes0.ndof)
                pdofs[:] = 0
                # FacetFESpace inner DOFs
                inner = prolM.GetInnerDofs(level)
                for j in range(dim):
                    pdofs[j * M.ndof:(j + 1) * M.ndof] = inner
                # pdofs[dim * M.ndof : dim * M.ndof + Q0.ndof] = 1
                # difficult to directly get global Q_l^T DOFs
                # he_prol, first perform he on M dofs, then p = Mass_p_inv * div(CR_M)
                pTransf = c_div * pFesEmb @ pMass_inv @ p_M_mix.mat.T @ MFesEmb.T
                pp.append([a0.mat.Inverse(pdofs, inverse="umfpack"), pTransf])
                # pp.append(a0.mat.Inverse(pdofs, inverse="umfpack"))
                # bk smoother
                # bjac = et.CreateSmoother(a0, {"blocktype": "vertexpatch"})
                pp.append(fes0.CreateSmoothBlocks(vertex=True, globalDofs=True))
                pp.append(bjac)
                pre.Update(a0.mat, pp)
        t2 = timeit.time()
        # estimate condition number
        if directSol:
            # lams = np.array([1, 1])
            inv = a0.mat.Inverse(fes0.FreeDofs(True), inverse='umfpack')
        else:
            # TODO: lams evaluation in mixed formulation?
            if uzawa:
                lams = EigenValues_Preconditioner(mat=a0.mat, pre=pre)
                if precond:
                    inv = CGSolver(a0.mat, pre, printrates=False, tol=1e-8, maxiter=100)
                else:
                    inv = IterSolver(mat=a0.mat, pre=pre, printrates=False, tol=1e-8, maxiter=200)
            else:
                # lams = np.array([1, 1])
                inv = GMResSolver(a0.mat, pre, printrates=False, tol=1e-10, maxiter=500)
                # inv = GMResSolver(a0.mat, a0.mat.Inverse(fes0.FreeDofs(True)), printrates=False, tol=1e-10, maxiter=500)
        t21 = timeit.time()
        # dirichlet BC
        # homogeneous Dirichlet assumed
        f.vec.data += a.harmonic_extension_trans * f.vec

        f0.vec.data += a0.harmonic_extension_trans * f0.vec
        it = 0
        if uzawa:
            # Augmented Lagrangian with Uzawa solver
            p_prev = GridFunction(Q)  # previous pressure solution
            p_prev.vec[:] = 0  # initialize
            for _ in range(uzIt):
                # dirichlet BC, TODO: Dirichlet BC set after each iteration?
                gfu0.vec.data[:] = 0
                gfu0.vec.data += inv * (f0.vec - MFesEmb @ p_M_mix.mat * p_prev.vec)
                p_prev.vec.data += c_div * (pMass_inv @ p_M_mix.mat.T @ MFesEmb.T * gfu0.vec.data)
                it = it + 1 if directSol else it + inv.iterations
            it /= uzIt
            # p_prev -> gfu
            gfu.vec.data += pFesEmb * p_prev.vec
        else:
            gfu0.vec.data += inv * f0.vec
            it = 1 if directSol else inv.iterations

        gfu.vec.data += Projector(fes0.FreeDofs(True), True) * gfu0.vec
        gfu.vec.data += a.harmonic_extension * gfu.vec
        gfu.vec.data += a.inner_solve * f.vec
        t3 = timeit.time()

        print(f"==> Time to find EIG Val: {t21 - t2:.1E}")
        print(f"==> COND: {max(lams)/min(lams):.2e}, MAX PREC LAM: {max(lams):.1E}; MIN PREC LAM: {min(lams):.1E}")
        print("==> Avg IT: %2.0i" % (it), "NDOFS: %.2e" % (
            sum(M.FreeDofs(True))*mesh.dim),
              "Assemble: %.2e Prec: %.2e Solve: %.2e" % (t1 - t0, t2 - t1, t3 - t21))
        # import netgen.gui
        # Draw(uh0, mesh)
        # input('continue?')

SolveBVP(0)
print('=====================')
print(
    f"DIM: {dim}, c_low: {c_lo:.1E}")
print('=====================')
print(f'level: {0}')
L2_uErr = sqrt(Integrate((uh - u_exact) * (uh - u_exact), mesh))
L2_LErr = sqrt(Integrate(InnerProduct((Lh - L_exact), (Lh - L_exact)), mesh))
L2_divErr = sqrt(Integrate(div(uh) * div(uh), mesh))
print(f"uh L2-error: {L2_uErr:.3E}")
print(f"Lh L2-error: {L2_LErr:.3E}")
print(f'uh divErr: {L2_divErr:.3E}, 1/epsilon: {c_div:.1E}')
print(f'uzawa solver: {bool(uzawa)}, uzIt: {uzIt}')
print(f'Direct Solver: {bool(directSol)}')
if not directSol: print(f"vertex-patch-GS steps: {ns}, var-V: {var}, W-cycle: {wc}")
print('==============================')
prev_uErr = L2_uErr
prev_uDivErr = L2_divErr
prev_LErr = L2_LErr
u_rate, uDiv_rate, L_rate = 0, 0, 0
level = 1
while True:
    with TaskManager():
        if dim == 2:
            mesh.ngmesh.Refine()
        elif dim == 3:
            mesh.Refine(onlyonce=True)
        meshRate = 2 if mesh.dim == 2 else sqrt(2)
        # mesh.ngmesh.Refine()
        # meshRate = 2
        # exit if total global dofs exceed a tol
        M.Update()
        if (M.ndof * dim > maxdofs):
            print(M.ndof * dim)
            break
        print(f'level: {level}')
        SolveBVP(level)
        # ===== convergence check =====
        L2_uErr = sqrt(Integrate((uh - u_exact) * (uh - u_exact), mesh))
        L2_LErr = sqrt(Integrate(InnerProduct((Lh - L_exact), (Lh - L_exact)), mesh))
        L2_divErr = sqrt(Integrate(div(uh) * div(uh), mesh))
        u_rate = log(prev_uErr / L2_uErr) / log(meshRate)
        L_rate = log(prev_LErr / L2_LErr) / log(meshRate)
        uDiv_rate = log(prev_uDivErr / L2_divErr) / log(meshRate)
        print(f"uh L2-error: {L2_uErr:.3E}, uh conv rate: {u_rate:.2E}")
        print(f"Lh L2-error: {L2_LErr:.3E}, Lh conv rate: {L_rate:.2E}")
        print(f'uh divErr: {L2_divErr:.3E}, uh div conv rate :{uDiv_rate:.2E}, 1/epsilon: {c_div:.1E}')
        print(f'uzawa solver: {bool(uzawa)}, uzIt: {uzIt}')
        print(f'Direct Solver: {bool(directSol)}')
        if not directSol: print(f"vertex-patch-GS steps: {ns}, var-V: {var}, W-cycle: {wc}")
        print('==============================')
        prev_uErr = L2_uErr
        prev_uDivErr = L2_divErr
        prev_LErr = L2_LErr
        level += 1
