from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *

from meshes import MakeTrig2DMesh

mesh = Mesh(unit_cube.GenerateMesh(maxh=2))

# fes = FESpace("BDM1", mesh)
fes = HDiv(mesh, order=2, loworderp1=True)
u,v = fes.TnT()

gfu = GridFunction(fes, nested=True)
# gfu.Set( (y,-x) )
gfu.Set( (y,x,x) )

# print ("coarse vec", gfu.vec)

print ("coarse vertices", mesh.nv)

Draw (gfu)
# input ("key")

a = BilinearForm(u*v*dx+div(u)*div(v)*dx)
c = Preconditioner(a, "multigrid", smoother="block", test=True, 
        inverse="umfpack")
a.Assemble()

#SetNumThreads(4)
with TaskManager(pajetrace=10**8):
  for l in range(3):
    mesh.Refine(onlyonce=False)
    print ("ndof = ", fes.ndof)
    gfu.Update()
    a.Assemble()

# hv = gfu.vec.CreateVector()
# hv2 = gfu.vec.CreateVector()
# hv.SetRandom()
# hv2.data = c * hv
# print (Norm(hv2))

# for face in mesh.faces:
#    f = face.nr
#    print ("face", f, " verts ", face.vertices, " ct = ", fes.couplingtype[3*f], #fes.couplingtype[3*f+1],fes.couplingtype[3*f+2],
#               hv2[3*f]) # , hv2[3*f+1], hv2[3*f+2])

# for f in range(mesh.nface):
#    print ("face", f, " has parents ", mesh.GetParentFaces(f))
    


# print (fes.loembedding)
# print (a.mat)
# print (a.loform.mat)

# am = a.mat
# alm = a.loform.mat
# em = fes.loembedding

# proj = em.CreateTranspose() @ am @ em
# print (proj)



