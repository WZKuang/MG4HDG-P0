from netgen.meshing import *
from netgen.csg import *
import ngsolve

def Make1DMesh(n, mapping = None, periodic=False):
    """
    Generate an equidistant 1D mesh with N cells

    Parameters
    ----------
    n : int
      Number of cells.

    mapping: lamda
      Mapping to transform the generated points. If None, the identity mapping is used.

    periodic: bool
      If True, the endpoints are identified to generate a periodic mesh.

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 1D NGSolve mesh

    """
    mesh = Mesh(dim=1)
    pids = []
    for i in range(n+1):
        x = i/n
        if mapping:
            x = mapping(x)
        pids.append (mesh.Add (MeshPoint(Pnt(x, 0, 0))))

    idx_inner = mesh.AddRegion("dom", dim=1)
    idx_left = mesh.AddRegion("left", dim=0)
    idx_right = mesh.AddRegion("right", dim=0)
        
    for i in range(n):
        mesh.Add(Element1D([pids[i],pids[i+1]],index=idx_inner))
    mesh.Add (Element0D( pids[0], index=idx_left))
    mesh.Add (Element0D( pids[n], index=idx_right))
    if periodic:
        mesh.AddPointIdentification(pids[0],pids[n],1,2)
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh

def MakeStructured2DMesh(quads=True, nx=10, ny=10, secondorder=False, periodic_x=False, periodic_y=False, mapping = None, bbpts=None, bbnames=None, flip_triangles=False):
    """
    Generate a structured 2D mesh

    Parameters
    ----------
    quads : bool
      If True, a quadrilateral mesh is generated. If False, the quads are split to triangles.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    secondorder : bool
      If True, second order curved elements are used.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    mapping: lamda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    bbpts : list
      List of points which should be handled as BBND and are named with bbnames. The mesh (nx, ny and mapping) must be constructed in such a way that the bbpts coincide with generated points. Otherwise an Exception is thrown.

    bbnames : list
      List of bbnd names as strings. Size must coincide with size of bbpts. Otherwise an Exception is thrown.

    flip_triangles : bool
      If set tot True together with quads=False the quads are cut the other way round

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 2D NGSolve mesh

    """
    mesh = Mesh()
    mesh.dim=2

    if (bbpts and bbnames) and len(bbpts) != len(bbnames):
        raise Exception("Lenght of bbnames does not coincide with length of bbpts!")

    found = []
    indbbpts = []
    if bbpts:
        for i in range(len(bbpts)):
            found.append(False)
            indbbpts.append(None)

    pids = []
    if periodic_y:
        minioni = []
        masteri = []
    if periodic_x:        
        minionj = []
        masterj = []
    for i in range(ny+1):
        for j in range(nx+1):
            x,y = j/nx, i/ny
            # if mapping:
            #    x,y = mapping(x,y)
            pids.append(mesh.Add (MeshPoint(Pnt(x,y,0))))
            if periodic_y:
                if i == 0:
                    minioni.append(pids[-1])
                if i == ny:
                    masteri.append(pids[-1])  
            if periodic_x:                       
                if j == 0:
                    minionj.append(pids[-1])
                if j == nx:
                    masterj.append(pids[-1])        
    if periodic_y:
        for i in range(len(minioni)):   
            mesh.AddPointIdentification(masteri[i],minioni[i],identnr=1,type=2)
    if periodic_x:            
        for j in range(len(minionj)):        
            mesh.AddPointIdentification(masterj[j],minionj[j],identnr=2,type=2)                                       

    # mesh.Add(FaceDescriptor(surfnr=1,domin=1,bc=1))
    idx_dom = mesh.AddRegion("dom", dim=2)
    idx_bottom = mesh.AddRegion("bottom", dim=1)
    idx_right  = mesh.AddRegion("right", dim=1)
    idx_top    = mesh.AddRegion("top", dim=1)
    idx_left   = mesh.AddRegion("left", dim=1)
    
    for i in range(ny):
        for j in range(nx):
            base = i * (nx+1) + j
            if quads:
                pnum = [base,base+1,base+nx+2,base+nx+1]
                elpids = [pids[p] for p in pnum]
                el = Element2D(idx_dom,elpids)
                if not mapping:
                    el.curved=False
                mesh.Add(el)
            else:
                if flip_triangles:
                    pnum1 = [base,base+1,base+nx+2]
                    pnum2 = [base,base+nx+2,base+nx+1]
                else:
                    pnum1 = [base,base+1,base+nx+1]
                    pnum2 = [base+1,base+nx+2,base+nx+1]
                elpids1 = [pids[p] for p in pnum1]
                elpids2 = [pids[p] for p in pnum2]
                mesh.Add(Element2D(idx_dom,elpids1)) 
                mesh.Add(Element2D(idx_dom,elpids2))                          

    for i in range(nx):
        mesh.Add(Element1D([pids[i], pids[i+1]], index=idx_bottom))
    for i in range(ny):
        mesh.Add(Element1D([pids[i*(nx+1)+nx], pids[(i+1)*(nx+1)+nx]], index=idx_right))
    for i in range(nx):
        mesh.Add(Element1D([pids[ny*(nx+1)+i+1], pids[ny*(nx+1)+i]], index=idx_top))
    for i in range(ny):
        mesh.Add(Element1D([pids[(i+1)*(nx+1)], pids[i*(nx+1)]], index=idx_left))

    # mesh.SetBCName(0, "bottom")        
    # mesh.SetBCName(1, "right")        
    # mesh.SetBCName(2, "top")        
    # mesh.SetBCName(3, "left")  

    mesh.Compress()       
    
    if secondorder:
        mesh.SecondOrder()
    
    if mapping:
        for p in mesh.Points():
            x,y,z = p.p
            x,y = mapping(x,y)
            p[0] = x
            p[1] = y

    for k in range(len(found)):
        i = 0
        for p in mesh.Points():
            if abs(p.p[0]-bbpts[k][0])+abs(p.p[1]-bbpts[k][1]) < 1e-6:
                indbbpts[k] = pids[i]
                found[k] = True
            i += 1
    for k in range(len(found)):
        if found[k] == False:
            raise Exception("bbpnt[",k,"] not in structured mesh!")

    for i in range(len(indbbpts)):
        mesh.Add(Element0D(indbbpts[i], index=i+1))
        mesh.SetCD2Name(i+1, bbnames[i])
    mesh.EnableTable ("parentedges")
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh

def MakeQuadMesh(nx=10, ny=10, periodic_x=False, periodic_y=False, mapping = None):
    """
    Generate a structured quadrilateral 2D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 2D NGSolve mesh

    """
    return MakeStructured2DMesh(quads=True, nx=nx, ny=ny, periodic_x=periodic_x, periodic_y=periodic_y, mapping=mapping)    

def MakeStructured3DMesh(hexes=True, nx=10, ny=None, nz=None, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False, prism=False):
    """
    Generate a structured quadrilateral 2D mesh

    Parameters
    ----------
    hexes: bool
      If True, a mesh consisting of hexahedra is generated.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    secondorder : bool
      If True, second order curved elements are used.
 
    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.

    prism : bool
      If True, a mesh consisting of prism is generated. If hexes and prism is set to True, also a mesh consisting of prism is generated.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh

    """
    if nz == None:
        if ny == None:
            nz = nx
        else:
            raise Exception("MakeStructured3DMesh: No default value for nz if nx and ny are provided")
    if ny == None:
        ny = nx
        
    netmesh = Mesh()
    netmesh.dim = 3

    if cuboid_mapping:
        P1 = (0,0,0)
        P2 = (1,1,1)
        if mapping:
            P1 = mapping(*P1)
            P2 = mapping(*P2)
        cube = OrthoBrick(Pnt(P1[0], P1[1], P1[2]), Pnt(P2[0], P2[1], P2[2])).bc(1)
        geom = CSGeometry()
        geom.Add(cube)
        netmesh.SetGeometry(geom)

    pids = []
    if periodic_x:
        minioni = []
        masteri = []
    if periodic_y:        
        minionj = []
        masterj = []
    if periodic_z:        
        minionk = []
        masterk = []        
    for i in range(nx+1):
        for j in range(ny+1):
            for k in range(nz+1):
                # x,y,z = mapping(i / nx, j / ny, k / nz)
                x,y,z = i / nx, j / ny, k / nz
                # if mapping:
                #   x,y,z = mapping(x,y,z)
                pids.append(netmesh.Add(MeshPoint(Pnt( x,y,z ))))
                if periodic_x:
                    if i == 0:
                        minioni.append(pids[-1])
                    if i == nx:
                        masteri.append(pids[-1])  
                if periodic_y:           
                    if j == 0:
                        minionj.append(pids[-1])
                    if j == ny:
                        masterj.append(pids[-1]) 
                if periodic_z:                    
                    if k == 0:
                        minionk.append(pids[-1])
                    if k == nz:
                        masterk.append(pids[-1])
    if periodic_x:
        for i in range(len(minioni)):   
            netmesh.AddPointIdentification(masteri[i],minioni[i],identnr=1,type=2)     
    if periodic_y:        
        for j in range(len(minionj)):            
            netmesh.AddPointIdentification(masterj[j],minionj[j],identnr=2,type=2) 
    if periodic_z:        
        for k in range(len(minionk)):            
            netmesh.AddPointIdentification(masterk[k],minionk[k],identnr=3,type=2)                                                      

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                base = i * (ny+1)*(nz+1) + j*(nz+1) + k
                baseup = base+(ny+1)*(nz+1)
                pnum = [base, base+1, base+(nz+1)+1, base+(nz+1),
                        baseup, baseup+1, baseup+(nz+1)+1, baseup+(nz+1)]
                if prism:
                    for qarr in [[0,3,4,1,2,5],[3,7,4,2,6,5]]:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element3D(1, elpids))
                elif hexes:
                    elpids = [pids[p] for p in pnum]
                    el = Element3D(1, elpids)
                    if not mapping:
                        el.curved = False
                    netmesh.Add(el)
                else:
                    #  a poor mans kuhn triangulation of a cube
                    for qarr in [[0, 4, 5, 6],
                                 [0, 6, 7, 4],
                                 [0, 3, 7, 6],
                                 [0, 1, 6, 5],
                                 [0, 1, 2, 6],
                                 [0, 3, 6, 2]]:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element3D(1, elpids))

    def AddSurfEls(p1, dxi, nxi, deta, neta, facenr):
        def add_seg(i, j, os):
            base = p1 + i*dxi + j*deta
            pnum = [base, base+os]
            elpids = [pids[p] for p in pnum]
            netmesh.Add(Element1D(elpids, index=facenr))
        for i in range(nxi):
            for j in [0,neta]:
                add_seg(i,j,dxi)
        for i in [0,nxi]:
            for j in range(neta):
                add_seg(i,j,deta)
        for i in range(nxi):
            for j in range(neta):
                base = p1 + i*dxi+j*deta
                pnum = [base, base+dxi, base+dxi+deta, base+deta]
                if prism:
                    if facenr <= 4:
                        qarr = [1,2,3,0]
                        elpids = [pids[pnum[p]] for p in qarr]
                        netmesh.Add(Element2D(facenr, elpids))
                    else:
                        qarrs = [[3, 1, 2], [3, 0, 1]]
                        for qarr in qarrs:
                            elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                            netmesh.Add(Element2D(facenr, elpids))
                elif hexes:
                    elpids = [pids[p] for p in pnum]
                    netmesh.Add(Element2D(facenr, elpids))
                else:
                    qarrs = [[0, 1, 2], [0, 2, 3]]
                    for qarr in qarrs:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element2D(facenr, elpids))

    #order is important!
    netmesh.Add(FaceDescriptor(surfnr=4, domin=1, bc=1))
    netmesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=2))
    netmesh.Add(FaceDescriptor(surfnr=5, domin=1, bc=3))
    netmesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=4))
    netmesh.Add(FaceDescriptor(surfnr=0, domin=1, bc=5))
    netmesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=6))
        
    # y-z-plane, smallest x-coord: ("back")
    AddSurfEls(0, 1, nz,  nz+1, ny, facenr=1) # y-z-plane
    # x-z-plane, smallest y-coord: ("left")
    AddSurfEls(0, (ny+1)*(nz+1), nx, 1, nz,facenr=2)
    # y-z-plane, largest x-coord: ("front")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -(nz+1), ny, -1, nz, facenr=3) 
    # x-z-plane, largest y-coord: ("right")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -1, nz, -(ny+1)*(nz+1), nx, facenr=4)
    # x-y-plane, smallest z-coord: ("bottom")
    AddSurfEls(0, nz+1, ny, (ny+1)*(nz+1), nx,facenr=5) 
    # x-y-plane, largest z-coord: ("top")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -(ny+1)*(nz+1), nx, -(nz+1), ny, facenr=6) 

    netmesh.SetBCName(0,"back")
    netmesh.SetBCName(1,"left")
    netmesh.SetBCName(2,"front")
    netmesh.SetBCName(3,"right")
    netmesh.SetBCName(4,"bottom")
    netmesh.SetBCName(5,"top")
    
    netmesh.Compress()

    if secondorder:
        netmesh.SecondOrder()
    
    if mapping:
        for p in netmesh.Points():
            x,y,z = p.p
            x,y,z = mapping(x,y,z)
            p[0] = x
            p[1] = y
            p[2] = z
            
    ngsmesh = ngsolve.Mesh(netmesh)
    # ngsmesh.ngmesh.Save("tmp.vol.gz")
    # ngsmesh = ngsolve.Mesh("tmp.vol.gz")
    return ngsmesh

def MakeHexMesh(nx=10, ny=10, nz=10, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False):
    """
    Generate a structured hexahedra 3D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh consisting of only hexahedra

    """
    return MakeStructured3DMesh(hexes=True, nx=nx, ny=ny, nz=nz, secondorder=secondorder, periodic_x=periodic_x, periodic_y=periodic_y, periodic_z=periodic_z, mapping=mapping, cuboid_mapping=cuboid_mapping)

def MakePrismMesh(nx=10, ny=None, nz=None, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False):
    """
    Generate a structured prism 3D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh consisting of only prism

    """
    return MakeStructured3DMesh(hexes=False, nx=nx, ny=ny, nz=nz, secondorder=secondorder, periodic_x=periodic_x, periodic_y=periodic_y, periodic_z=periodic_z, mapping=mapping, cuboid_mapping=cuboid_mapping, prism=True)

def MakeTet3DMesh(nref=0):
    netmesh = Mesh()
    netmesh.dim = 3

    pids = []
    pids.append(netmesh.Add(MeshPoint(Pnt( 1,0,0 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,1,0 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,0,1 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,0,0 ))))

    # one tet
    netmesh.Add(Element3D(1, [1,2,3,4]))
    
    # four faces
    netmesh.Add(FaceDescriptor(surfnr=0, domin=1, bc=1))
    netmesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=2))
    netmesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=3))
    netmesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=4))

    # diagonal-plane
    netmesh.Add(Element1D([2,3], index=1))
    netmesh.Add(Element1D([3,4], index=1))
    netmesh.Add(Element1D([4,2], index=1))
    netmesh.Add(Element2D(1, [2,3,4]))
    # x-plane
    netmesh.Add(Element1D([3,4], index=2))
    netmesh.Add(Element1D([4,1], index=2))
    netmesh.Add(Element1D([1,3], index=2))
    netmesh.Add(Element2D(2, [3,4,1]))
    # y-plane
    netmesh.Add(Element1D([4,1], index=3))
    netmesh.Add(Element1D([1,2], index=3))
    netmesh.Add(Element1D([2,4], index=3))
    netmesh.Add(Element2D(3, [4,1,2]))
    # z-plane
    netmesh.Add(Element1D([1,2], index=4))
    netmesh.Add(Element1D([2,3], index=4))
    netmesh.Add(Element1D([3,1], index=4))
    netmesh.Add(Element2D(4, [1,2,3]))

    netmesh.SetBCName(0,"diag")
    netmesh.SetBCName(1,"xplane")
    netmesh.SetBCName(2,"yplane")
    netmesh.SetBCName(3,"zplane")
    
    netmesh.Compress()

    ngsmesh = ngsolve.Mesh(netmesh)
    for i in range(nref):
        ngsmesh.Refine()
    return ngsmesh

def MakeTet3DMeshB(bisect=0):
    netmesh = Mesh()
    netmesh.dim = 3

    pids = []
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,0,0 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 1,0,0 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,1,0 ))))
    pids.append(netmesh.Add(MeshPoint(Pnt( 0,0,1 ))))
    # midpoint between 0 and 1
    pids.append(netmesh.Add(MeshPoint(Pnt( 0.5,0,0 ))))
    # two tets
    netmesh.Add(Element3D(1, [1,5,3,4]))
    netmesh.Add(Element3D(1, [2,3,4,5]))
    
    # four faces
    netmesh.Add(FaceDescriptor(surfnr=0, domin=1, bc=1))
    netmesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=2))
    netmesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=3))
    netmesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=4))

    # diagonal-plane
    netmesh.Add(Element1D([2,3], index=1))
    netmesh.Add(Element1D([3,4], index=1))
    netmesh.Add(Element1D([2,4], index=1))
    netmesh.Add(Element2D(1, [2,3,4]))
    # x-plane
    netmesh.Add(Element1D([3,4], index=2))
    netmesh.Add(Element1D([1,4], index=2))
    netmesh.Add(Element1D([1,3], index=2))
    netmesh.Add(Element2D(2, [1, 3,4]))
    # y-plane
    netmesh.Add(Element1D([1,4], index=3))
    netmesh.Add(Element1D([1,5], index=3))
    netmesh.Add(Element1D([2,5], index=3))
    netmesh.Add(Element1D([2,4], index=3))
    netmesh.Add(Element1D([4,5], index=3))
    netmesh.Add(Element2D(3, [1, 4, 5]))
    netmesh.Add(Element2D(3, [2, 4, 5]))
    # z-plane
    netmesh.Add(Element1D([1,5], index=4))
    netmesh.Add(Element1D([2,5], index=4))
    netmesh.Add(Element1D([2,3], index=4))
    netmesh.Add(Element1D([1,3], index=4))
    netmesh.Add(Element1D([3,5], index=4))
    netmesh.Add(Element2D(4, [1,3,5]))
    netmesh.Add(Element2D(4, [2,3,5]))

    netmesh.SetBCName(0,"diag")
    netmesh.SetBCName(1,"xplane")
    netmesh.SetBCName(2,"yplane")
    netmesh.SetBCName(3,"zplane")
    
    netmesh.Compress()

    ngsmesh = ngsolve.Mesh(netmesh)
    for i in range(nref):
        ngsmesh.Refine()
    return ngsmesh

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

def MakeTrig2DMeshR(nref=0, y0=1):
    mesh = Mesh()
    mesh.dim=2
    pids=[]
    pids.append(mesh.Add (MeshPoint(Pnt(0,y0,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(0,0,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(1,0,0))))
              
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

def MakeTrig2DMeshX():
    mesh = Mesh()
    mesh.dim=2
    pids=[]
    pids.append(mesh.Add (MeshPoint(Pnt(1,0,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(0,1,0))))
    pids.append(mesh.Add (MeshPoint(Pnt(0,0,0))))
              
    # one trig    
    idx_dom = mesh.AddRegion("dom", dim=2)
    mesh.Add(Element2D(idx_dom, [1,2,3]))
    idx_bottom = mesh.AddRegion("bottom", dim=1)
    mesh.Add(Element1D([1,2], index=idx_bottom))
    mesh.Add(Element1D([2,3], index=idx_bottom))
    mesh.Add(Element1D([3,1], index=idx_bottom))
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh
