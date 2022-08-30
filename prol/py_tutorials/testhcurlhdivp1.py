from netgen.csg import unit_cube
from ngsolve import *
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

hdiv = True

if hdiv:
    fes = FESpace("BDM1", mesh)
else:
    fes = FESpace("HCurlP1", mesh)
u,v = fes.TnT()

gfu = GridFunction(fes, nested=True)
if hdiv:
    gfu.Set( (z,x, y) )
else:
    gfu.Set( (x, y, z) )

if hdiv:
    a = BilinearForm(u*v*dx+div(u)*div(v)*dx)
else:
    a = BilinearForm(u*v*dx+curl(u)*curl(v)*dx)
c = Preconditioner(a, "multigrid", smoother="block", 
        updateall=False, 
        coarsetype="direct",
        coarsesmoothingsteps=1, # ??? 
        updatealways=True,
        test=True)
a.Assemble()

for l in range(5):
    mesh.Refine()
    gfu.Update()
    a.Assemble()
    print ("ndof = ", fes.ndof)


# Draw (gfu[0], mesh, "ux")

