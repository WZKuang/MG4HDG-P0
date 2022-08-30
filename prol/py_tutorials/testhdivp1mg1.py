from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from ngsolve.krylovspace import CGSolver
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

fes = FESpace("BDM1", mesh)
u,v = fes.TnT()

gfu = GridFunction(fes, nested=True)
# gfu.Set( (y,-x) )
gfu.Set( (y,x,x) )

# print ("coarse vec", gfu.vec)

print ("coarse vertices", mesh.nv)

Draw (gfu)
# input ("key")

a = BilinearForm(u*v*dx+div(u)*div(v)*dx)
c = Preconditioner(a, "multigrid", smoother="block", test=True)
a.Assemble()

f = LinearForm(v[0]*dx)
with TaskManager():
  for l in range(4):

    mesh.Refine() # onlyonce=True)
    gfu.Update()
    a.Assemble()
    print ("ndof = ", fes.ndof)
    #f.Assemble()
    #inv = CGSolver(a.mat, c.mat, printing=True)
    #gfu.vec.data = inv*f.vec


stop

hv = gfu.vec.CreateVector()
hv2 = gfu.vec.CreateVector()
hv.SetRandom()
hv2.data = c * hv

# print (hv2)

for face in mesh.faces:
    f = face.nr
    print ("face", f, " verts ", face.vertices, " ct = ", fes.couplingtype[3*f], #fes.couplingtype[3*f+1],fes.couplingtype[3*f+2],
               hv2[3*f]) # , hv2[3*f+1], hv2[3*f+2])

for f in range(mesh.nface):
    print ("face", f, " has parents ", mesh.GetParentFaces(f))
    
