from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from ngsolve.meshes import *
from meshes import *
import netgen.gui
#from prol import *
#mesh = MakeStructured2DMesh(quads=False, nx=1, ny=1)
mesh = MakeTrig2DMeshX()
#mesh = MakeTet3DMesh()
i = ElementId(VOL,0)
trafo = mesh.GetTrafo(i)
ir = IntegrationRule(TRIG, 2)
print(ir)
for p in ir:
    print("ref: (", p.point[0],",",
            p.point[1],",", p.point[2], ")")
    print("phy: (", trafo(p).point[0],",",
            trafo(p).point[1],")")
Draw(mesh)
input("?")
