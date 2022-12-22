from ngsolve import *
from netgen.geom2d import SplineGeometry
import matplotlib.pyplot as plt
from prol import *
from mymg import *
from ngsolve.la import EigenValues_Preconditioner
from ngsolve.krylovspace import CGSolver


# ### Geometry
# 
# The following geometry represents a heated chip embedded in another material that conducts away the heat.


#   point numbers 0, 1, ... 11
#   sub-domain numbers (1), (2), (3)
#  
#
#             7-------------6
#             |             |
#             |     (2)     |
#             |             |
#      3------4-------------5------2
#      |                           |
#      |             11            |
#      |           /   \           |
#      |         10 (3) 9          |
#      |           \   /     (1)   |
#      |             8             |
#      |                           |
#      0---------------------------1
#

def MakeGeometry():
    geometry = SplineGeometry()
    
    # point coordinates ...
    pnts = [ (0,0), (1,0), (1,0.6), (0,0.6),              (0.2,0.6), (0.8,0.6), (0.8,0.8), (0.2,0.8),              (0.5,0.15), (0.65,0.3), (0.5,0.45), (0.35,0.3) ]
    pnums = [geometry.AppendPoint(*p) for p in pnts]
    
    # start-point, end-point, boundary-condition, left-domain, right-domain:
    lines = [ (0,1,1,1,0), (1,2,2,1,0), (2,5,2,1,0), (5,4,2,1,2), (4,3,2,1,0), (3,0,2,1,0),               (5,6,2,2,0), (6,7,2,2,0), (7,4,2,2,0),               (8,9,2,3,1), (9,10,2,3,1), (10,11,2,3,1), (11,8,2,3,1) ]
        
    for p1,p2,bc,left,right in lines:
        geometry.Append(["line", pnums[p1], pnums[p2]], bc=bc, leftdomain=left, rightdomain=right)

    geometry.SetMaterial(1,"base")
    geometry.SetMaterial(2,"top")
    geometry.SetMaterial(3,"chip")    

    return geometry

from netgen.csg import *
def MakeGeometry3D():
    geo = CSGeometry()
    box = OrthoBrick(Pnt(0,0,0), Pnt(1,0.6,0.6)).mat("base")
    left  = Plane (Pnt(0.5,0.15,0), Vec(-1,-1,0) )
    right = Plane (Pnt(0.5,0.45,1), Vec( 1,1,0) )
    front = Plane (Pnt(0.5,0.15,0), Vec(1,-1,0) )
    back  = Plane (Pnt(0.5,0.45,1), Vec(-1, 1,0) )
    bot   = Plane (Pnt(0,0,0.225), Vec(0,0,-1) )
    top   = Plane (Pnt(1,1,0.375), Vec(0,0, 1) )
    box2 = (left*right*front*back*bot*top).mat("chip")

#     box2 = Sphere(Pnt(0.5,0.3,0.3), 0.1).mat("chip")
#     box2 = OrthoBrick(Pnt(0.4,0.2,0.2), Pnt(0.6,0.4,0.3)).mat("chip")
    box3 = OrthoBrick(Pnt(0.2,0.6,0.2), Pnt(0.8,0.8,0.4)).mat("top")
    geo.Add(box-box2)
    geo.Add(box3)
    geo.Add(box2)
    return geo


# ==========coefficients==========
import sys
if len(sys.argv) < 5:
    print("NOT ENOUGH ARGUMENTS!!! dim + sm steps + block? + low order coef")
    exit(1)
dim = int(sys.argv[1])
smStep = int(sys.argv[2])
block = bool(int(sys.argv[3]))
rac_cof = int(sys.argv[4])
drawResults = True
adaptive = True
# ===============================
maxLevel = 25
if dim == 2:
    mesh = Mesh(MakeGeometry().GenerateMesh(maxh=1/4))
    index = 1
    maxdofs=1e6
else:
    mesh = Mesh(MakeGeometry3D().GenerateMesh(maxh=1/4))
    index = 12
    maxdofs=5e6
ne = mesh.ne

# HDG
V = VectorL2(mesh) 
W = L2(mesh, order=1)
M = FacetFESpace(mesh, dirichlet=[index])
fes = M*V*W
et = meshTopology(mesh, mesh.dim)
et.Update() 
if mesh.dim==2:
    prol = FacetProlongationTrig(mesh, et) 
else:
    prol = FacetProlongationTet(mesh, et) 

n = specialcf.normal(mesh.dim) 
h = specialcf.mesh_size
a = BilinearForm(fes, symmetric=False, condense=True)
# non-variant coefficient auxiliary operator
a_ax = BilinearForm(fes, symmetric=False, condense=True)
(uhat, q, u),(vhat,r, v) = fes.TnT()  

# one heat conductivity coefficient per sub-domain
lam = CoefficientFunction([1, 1000, 10])
alpha = 1
c_rac = CoefficientFunction([rac_cof, rac_cof, rac_cof])

if mesh.dim == 3:
    ir_c = IntegrationRule(points = [(1/3, 1/3, 0), (1/3, 0, 1/3), (0, 1/3, 1/3), (1/3, 1/3, 1/3)],
                           weights = [1/16, 1/16, 1/16, 1/16])
    a += (1/lam * q * r + c_rac * u * v) * dx(intrules = {QUAD:ir_c})
    a_ax += (1/lam * q * r) * dx(intrules = {QUAD:ir_c})
elif mesh.dim == 2:
    ir_c = IntegrationRule(points = [(0.5, 0), (0, 0.5), (0.5, 0.5)],
                           weights= [1/6, 1/6, 1/6])
    a += (1/lam * q * r + c_rac * u * v) * dx(intrules = {TRIG:ir_c})
    a_ax += (1/lam * q * r) * dx(intrules = {TRIG:ir_c})

if mesh.dim ==2:
    ir = IntegrationRule(SEGM, 1) 
    a += (-uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                          intrules={SEGM:ir})
    a_ax += (c_rac * h / 3 * uhat * vhat
             -uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                          intrules={SEGM:ir})
else:
    ir = IntegrationRule(TRIG, 1) 
    a += (-uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                          intrules={TRIG:ir})
    a_ax += (c_rac * h / 4 * uhat * vhat
             -uhat*r*n+q*n*vhat
          + lam*alpha/h*(u-uhat)*(v-vhat))*dx(element_boundary=True,
                                          intrules={TRIG:ir})
# heat-source in inner subdomain
f = LinearForm(fes)
f += CoefficientFunction([0, 0, 1])*v * dx
a.Assemble()
# a_ax.Assemble()
# =========preconditioner init==============
pre = MultiGrid(a.mat, prol, nc=M.ndof,
                coarsedofs=fes.FreeDofs(True), w1=0.6, 
                nsmooth=smStep, sm="gs", var=False)

gfu = GridFunction(fes)
uhath, qh, uh = gfu.components
# vtk = VTKOutput(ma=mesh,coefs=[uh], names=["sol"], filename="adpt"+str(dim),subdivision=0)

def SolveBVP():
    with TaskManager():
        fes.Update()
        gfu.Update()
        a.Assemble()
        # a_ax.Assemble()
        f.Assemble()
        if mesh.ne > ne:
            et.Update()
            # ===========point or block============
            # block
            if block:
                pp = [fes.FreeDofs(True), M.ndof, [], VertexPatchBlocks(mesh,M)]
            # point
            else:
                pp = [fes.FreeDofs(True), M.ndof, [], fes.FreeDofs(True)]
            pre.Update(a.mat, pp)
        
        lams = EigenValues_Preconditioner(mat=a.mat, pre=pre)
        inv = CGSolver(a.mat, pre, printing=False, tol=1e-8, maxiter=80)
    #     inv = a.mat.Inverse(fes.FreeDofs(True))
        f.vec.data += a.harmonic_extension_trans * f.vec  
        gfu.vec.data = inv*f.vec
        it = inv.iterations   
        gfu.vec.data += a.harmonic_extension * gfu.vec 
        gfu.vec.data += a.inner_solve * f.vec
        print(f"IT: {it}, cond: {max(lams)/min(lams):.2e}, Facet DOFS: {M.ndof}")
        
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
    print(f"ndof = {fes.ndof}, global dofs: {sum(fes.FreeDofs(True))}  maxerr = {maxerr**0.5:.2e}")
    
    # mark for refinement:
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, eta2[el.nr] > 0.25 * maxerr)


print('==============================')
print(f'DIM: {dim}, Smooth Steps: {smStep}, Block GS: {block}, Reaction Coef: {rac_cof}')
print('==============================')
level = 0
SolveBVP()
while M.ndof < maxdofs and level < maxLevel:  
    level += 1
    print('===============')
    print(f'Level: {level}')
    if adaptive:  CalcError()
    # mesh.Refine()
    if dim == 2:
        if adaptive:
            mesh.Refine()
        else:
            mesh.ngmesh.Refine()
    else:
        mesh.Refine(onlyonce=True)
        # if adaptive:
        #     mesh.Refine()
        # else:
        #     mesh.Refine(onlyonce = True)
    SolveBVP()

