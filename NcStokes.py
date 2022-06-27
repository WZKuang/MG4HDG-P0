from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import *
from mymg import *
from ngsolve.la import EigenValues_Preconditioner
# ============ parameters ===========
c_vis, c_div = 1, 1e1
c_lo = 10  # int(sys.argv[3])
dim = 2

# direct = True
# ==================================

if dim == 2:
    mesh = Mesh(unit_square.GenerateMesh(maxh=1 / 4))
    maxdofs = 1e5
else:
    mesh = Mesh(unit_cube.GenerateMesh(maxh=1 / 3))
    maxdofs = 1e5

# ========= NC space
V = FESpace('nonconforming', mesh, dirichlet='.*')
Q = L2(mesh, lowest_order_wb = False)
# V = H1(mesh, order = 2, dirichlet = '.*')
# Q = H1(mesh, order = 1)#, lowest_order_wb=mixed)


# ======= AL-uzawa solver
fes = V * V
gfu = GridFunction(fes)
a = BilinearForm(fes)
ir_c = IntegrationRule(points=[(0.5, 0), (0, 0.5), (0.5, 0.5)],
                       weights=[1 / 6, 1 / 6, 1 / 6])
(ux, uy), (vx, vy) = fes.TnT()
u = CF((ux, uy))
v = CF((vx, vy))
GradU = CF((grad(ux), grad(uy)))
GradV = CF((grad(vx), grad(vy)))
divU = grad(ux)[0] + grad(uy)[1]
divV = grad(vx)[0] + grad(vy)[1]
a += (c_vis * InnerProduct(GradU, GradV) + c_lo * u * v + c_div * divU * divV) * dx(intrules={TRIG: ir_c})
uhx, uhy = gfu.components
uh = CF((uhx, uhy))
# ===== Uzawa operator B
pNC, qNC = Q.TnT()
# grad mat: p -> NC * dim
p_NC_mix = BilinearForm(trialspace=Q, testspace=fes)
p_NC_mix += pNC * divV * dx
pMass = BilinearForm(Q)
pMass += pNC * qNC * dx


# ========= Saddle point problem
fes0 = V * V * Q
a0 = BilinearForm(fes0)
gfu0 = GridFunction(fes0)
(ux0, uy0, p0), (vx0, vy0, q0) = fes0.TnT()
u0 = CF((ux0, uy0))
v0 = CF((vx0, vy0))
GradU0 = CF((grad(ux0), grad(uy0)))
GradV0 = CF((grad(vx0), grad(vy0)))
divU0 = grad(ux0)[0] + grad(uy0)[1]
divV0 = grad(vx0)[0] + grad(vy0)[1]
a0 += (c_vis * InnerProduct(GradU0, GradV0) + c_lo * u0 * v0 + p0 * divV0 + q0 * divU0) * dx(intrules={TRIG: ir_c})
uhx0, uhy0, ph0 = gfu0.components
uh0 = CF((uhx0, uhy0))

# ===== exact solution
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

p_exact = -100 * x * (1 - x) * (1 - y) - 100 / 12


# ===== linear forms
f = LinearForm(fes)
f += (-c_vis * (4 * y * (1 - y) * (2 * y - 1) * ((1 - 2 * x) ** 2 - 2 * x * (1 - x))
                + 12 * x ** 2 * (1 - x) ** 2 * (1 - 2 * y))
      + 100 *(1 - 2 * x) * (1 - y)) * v[0] * dx(intrules={TRIG: ir_c})
f += (-c_vis * (4 * x * (1 - x) * (1 - 2 * x) * ((1 - 2 * y) ** 2 - 2 * y * (1 - y))
                + 12 * y ** 2 * (1 - y) ** 2 * (2 * x - 1))
      - 100 * x * (1 - x)) * v[1] * dx(intrules={TRIG: ir_c})
f += c_lo * u_exact * v * dx(intrules={TRIG: ir_c})

f0 = LinearForm(fes0)
f0 += (-c_vis * (4 * y * (1 - y) * (2 * y - 1) * ((1 - 2 * x) ** 2 - 2 * x * (1 - x))
                + 12 * x ** 2 * (1 - x) ** 2 * (1 - 2 * y))
      + 100 * (1 - 2 * x) * (1 - y)) * v0[0] * dx(intrules={TRIG: ir_c})
f0 += (-c_vis * (4 * x * (1 - x) * (1 - 2 * x) * ((1 - 2 * y) ** 2 - 2 * y * (1 - y))
                + 12 * y ** 2 * (1 - y) ** 2 * (2 * x - 1))
      - 100 * x * (1 - x)) * v0[1] * dx(intrules={TRIG: ir_c})
f0 += c_lo * u_exact * v0 * dx(intrules={TRIG: ir_c})


def SolveBVP(level):
    fes.Update(); fes0.Update()
    gfu.Update(); gfu0.Update()
    a.Assemble(); a0.Assemble()
    f.Assemble(); f0.Assemble()

    Q.Update()
    p_NC_mix.Assemble()
    pMass.Assemble()
    pMass_inv = pMass.mat.CreateSmoother()
    Proj_bd = Projector(fes.FreeDofs(), True)

    # ===== p L2 norm
    Proj_p = Embedding(fes0.ndof, fes0.Range(2)).T  # fes0 -> pNC
    # dirichlet BC
    # homogeneous Dirichlet assumed
    # ===== True solution, incompressible
    gfu0.vec.data += a0.mat.Inverse(fes0.FreeDofs(), inverse = 'umfpack') * f0.vec
    # ===== AL-uzawa solver, epsilon = 1/c_div
    pTmp = GridFunction(Q)
    pTmp.vec.data[:] = 0
    errP = pTmp.vec.CreateVector()
    errP0 = pTmp.vec.CreateVector()
    errP0.data = Proj_p * gfu0.vec - pTmp.vec
    for _ in range(8):
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse='umfpack') * \
                        (f.vec - p_NC_mix.mat * pTmp.vec)
        pTmp.vec.data += c_div * pMass_inv @ p_NC_mix.mat.T * gfu.vec
        errP.data = Proj_p * gfu0.vec - pTmp.vec
        # print(f'err P norm: {Norm(errP):.2E}, Err_P norm decrease rate: {Norm(errP) / Norm(errP0):.1E}')
        errP0.data = errP
        # Draw(p_exact, mesh, 'p')
        # input('true p?')
        # Draw(ph0, mesh, 'p')
        # input('saddle p?')
        # Draw(pTmp, mesh, 'p')
        # input('AL-uzawa p?')

# SolveBVP(0)
u_rate = 0
level = 1
while True:
    with TaskManager():
        mesh.ngmesh.Refine()
        # exit if total global dofs exceed a tol
        V.Update()
        if (V.ndof * dim > maxdofs):
            print(V.ndof * dim)
            break
        if level == 1:
            print('=====================')
            print(
                f"DIM: {dim}, c_vis: {c_vis:.1E}, c_low: {c_lo:.1E}, 1/eps: {c_div:.1E}")
            print('=====================')
        print(f'level: {level}')
        SolveBVP(level)
        # ===== convergence check =====
        meshRate = 2
        L2_uErr = sqrt(Integrate((uh - u_exact) * (uh - u_exact), mesh))
        divUh = grad(uhx)[0] + grad(uhy)[1]
        L2_divErr = sqrt(Integrate(divUh * divUh, mesh))
        if level != 1:
            u_rate = log(prev_uErr / L2_uErr) / log(meshRate)
        print(f"uh L2-error: {L2_uErr:.3E}, uh conv rate: {u_rate:.2E}")
        print(f"uh divErr: {L2_divErr:.1E}, 1/epsilon: {c_div:.1E}")
        print('==============================')
        prev_uErr = L2_uErr
        level += 1

