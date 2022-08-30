from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.krylovspace import CGSolver
ngsglobals.msg_level=0

SetTestoutFile ("test1.out")

def checkMat(compound):
    mesh = unit_square.GenerateMesh(maxh=1)
    mesh.EnableTable ("parentedges")
    mesh = Mesh(mesh)
    if compound:
        V = FESpace("HCurlP1", mesh, dirichlet=".*", dgjumps=True)
        W = H1(mesh, order = 2, dirichlet=".*")
        fes = FESpace([V,W], dgjumps=True)
        (u, p), (v, q) = fes.TnT()
    else:
        fes = FESpace("HCurlP1", mesh, dirichlet=".*", dgjumps=True)
        u, v = fes.TnT()
    gfu = GridFunction(fes)
    a = BilinearForm(fes, printelmat=True)
    
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)
    
    # Normal velocity
    def norm(v):
        return (v*n)*n
    
    # viscous term (laplacian)
    # a += InnerProduct(Grad(u), Grad(v))*dx
    
    avgdu = 0.5*(Grad(u)+Grad(u.Other()))*n
    avgdv = 0.5*(Grad(v)+Grad(v.Other()))*n
    jmpu = norm(u-u.Other())
    jmpv = norm(v-v.Other())
    # a += (-avgdu*jmpv-avgdv*jmpu+10/h*jmpu*jmpv)*dx(skeleton=True)
    a += (10/h*jmpu*jmpv)*dx(skeleton=True)
    
    a.Assemble()
    print(a.mat)
    #for k in range(1):
    #    mesh.Refine()
    #    fes.Update()
    #    gfu.Update()
    #    a.Assemble()
    #    rows,cols,vals = a.mat.COO()
    #    # import scipy.sparse as sp
    #    # A = sp.csr_matrix((vals, (rows, cols)))
    #    # print(A[10,11])
    #    print (a.mat)

print("single HCurlP1:")
checkMat(False) 
print("HCurlP1 in a compound space:")
checkMat(True) 

