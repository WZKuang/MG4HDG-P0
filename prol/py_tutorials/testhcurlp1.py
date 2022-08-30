# test hdiv/hcurl MG
from netgen.meshing import *
import ngsolve
def MakeTrig2DMesh(nref=0, y0=1):
    mesh = Mesh()
    mesh.dim=2
    pids=[]
    pids.append(mesh.Add (MeshPoint(Pnt(0,0,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(1,0,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(0,y0,0))))
              
    # one trig    
    idx_dom = mesh.AddRegion("dom", dim=2)
    mesh.Add(Element2D(idx_dom, [1,2,3]))
    idx_bottom = mesh.AddRegion("bottom", dim=1)
    mesh.Add(Element1D([1,2], index=idx_bottom))
    mesh.Add(Element1D([2,3], index=idx_bottom))
    mesh.Add(Element1D([3,1], index=idx_bottom))
    ngsmesh = ngsolve.Mesh(mesh)
    for i in range(nref):
        ngsmesh.Refine()
    return ngsmesh

mesh = MakeTrig2DMesh()

from netgen.csg import unit_cube
from netgen.geom2d import unit_square
from ngsolve import *

hdiv = True
bisect = False
if hdiv:
    fes = FESpace("BDM1", mesh)
else:
    fes = FESpace("HCurlP1", mesh)

gfu = GridFunction(fes, nested=True)
if hdiv:
    gfu.Set( (y, x) )
else:
    gfu.Set( (x, y) )
    #gfu.Set( (y, -x) )

u,v = fes.TnT()
if hdiv:
    a = BilinearForm(u*v*dx+div(u)*div(v)*dx)
else:
    a = BilinearForm(u*v*dx+curl(u)*curl(v)*dx)
c = Preconditioner(a, "multigrid", smoother="block", test=True)
a.Assemble()

for l in range(7):
    if bisect:
        mesh.Refine()
    else:
        mesh.ngmesh.Refine()
    gfu.Update()
    a.Assemble()
    print ("ndof = ", fes.ndof)

Draw (gfu[0], mesh, "ux")

