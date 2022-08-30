#pragma once

/// from ngxfem
#include <comp.hpp>
#include <multigrid.hpp>

namespace ngmg
{

  // mesh preprocessing
  //
  // GetVertLevel (int level)::: # of total vertices up to level [level]
  // VertLevel (int vtnr): the mesh lvl for a vertex [vtnr]
  // GetEdgeLevel (int level)::: # of total edges up to level [level]
  // EdgeNodes (int ednr): vertex #s of an edge [ednr]

  // ************************EDGE DATA STRUCTURE*****************************
  //                 [THIS IS NEEDED FOR P2-H1/P1-HCurl prol]
  // GetEdgeDofNumber (int lvl): edge dofs ordering of ALL edges at mesh level [lvl]
  // [We first order all active edges, then append those inactive edges]
  //
  //
  //
  // ************************FACE DATA STRUCTURE*****************************
  //                 [THIS IS NEEDED FOR P1-HDiv tet prol]
  // GetFaceDofNumber (int lvl): face dofs ordering of ALL faces (including intermediate faces) at mesh level [lvl]
  // [We first order all active faces, then append those inactive faces]
  //
  class meshTopology // Trig/Tet only
  {
    shared_ptr<MeshAccess> ma;
    int dim; // dimension
    
    Array<int> nvlevel; // # of vertices up to level
    Array<int> nelevel; // # of edges up to  level
    Array<int> nflevel; // # of ACTIVE facets up to level
    Array<int> nclevel; // # of elements up to  level
    Array<short int> vertexlevel, vertexlevel0;
    Array<short int> edgelevel, edgelevel0;
    Array<short int> facelevel, facelevel0;
    Array<int> faceIC, edgeIC; // in active edges/faces
    Array<INT<12>> ccedges, ccedges0;
    Array<INT<2>> ceedges, cefaces; // coarse element #s
    Array<Array<int>> ccfaces;
    Array<Array<float>> cwfaces;
    Array<INT<2>> edgenodes;
    Array<INT<3>> facenodes;
    Array<shared_ptr<Array<int>>> edgeDofNumber, faceDofNumber;
    Array<Array<INT<3>>> elemEdges, elemVerts;
    Array<Array<INT<2>>> edgeElems, edgeElemsLvl;

    Array<Array<INT<4>>> elemFaces, elemVerts3D;
    Array<Array<INT<2>>> faceElems, faceElemsLvl;

    public:
    meshTopology(shared_ptr<MeshAccess> ama, int adim): ma(ama), dim(adim){ ; }
    virtual ~meshTopology (){;}
    void Update();
    
    virtual shared_ptr<BaseBlockJacobiPrecond> CreateSmoother (shared_ptr<BilinearForm> bfa, 
        Flags & flags) const;
    // TODO 
    virtual shared_ptr<BaseBlockJacobiPrecond> CreateInnerSolve (shared_ptr<BilinearForm> bfa, 
        Flags & flags, int finelevel) const;

    shared_ptr<MeshAccess> GetMeshAccess() const { return ma; }

    size_t GetVertLevel (int level) const {return nvlevel[level];}
    size_t GetEdgeLevel (int level) const {return nelevel[level];}
    size_t GetFaceLevel (int level) const {return nflevel[level];}
    size_t GetCellLevel (int level) const {return nclevel[level];}

    auto CCEdges (int ednr) const { return ccedges[ednr]; }
    auto CCFaces (int ednr, int j) const { return ccfaces[ednr][j]; }
    auto CWFaces (int ednr, int j) const { return cwfaces[ednr][j]; }

    auto EdgeNodes (int ednr) const { return edgenodes[ednr];}
    auto VertLevel (int vtnr) const { return vertexlevel[vtnr];}
    auto GetEdgeDofNumber (int lvl) const { return edgeDofNumber[lvl];}
    auto GetFaceDofNumber (int lvl) const { return faceDofNumber[lvl];}
    auto GetEdgeIC (int lvl) const { return edgeIC[lvl];}
    auto GetFaceIC (int lvl) const { return faceIC[lvl];}
  };

  // Facet Prolongation for FacetFESpace
  class FacetProlongationTrig : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const meshTopology& et;
    public:
    FacetProlongationTrig(shared_ptr<MeshAccess> ama, const meshTopology & aet)
      : ma(ama), et(aet){;} 
    virtual ~FacetProlongationTrig() { ; }

    virtual void Update (const FESpace & fes){;} // no need to update 
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("FacetProlongationTrig::CreateProlongationMatrix not implemented!");
    }
    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const;
    
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };
  
  // Facet Prolongation for nonconforming FEM
  class FacetProlongationTrig2 : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const meshTopology& et;
    public:
    FacetProlongationTrig2(shared_ptr<MeshAccess> ama, const meshTopology & aet)
      : ma(ama), et(aet){;} 
    virtual ~FacetProlongationTrig2() { ; }

    virtual void Update (const FESpace & fes){;} // no need to update 
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("FacetProlongationTrig::CreateProlongationMatrix not implemented!");
    }
    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const;
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };
  
  
  // Facet Prolongation for FacetFESpace
  class FacetProlongationTet : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const meshTopology& et;
    public:
    FacetProlongationTet(shared_ptr<MeshAccess> ama, const meshTopology & aet)
      : ma(ama), et(aet){;} 
    virtual ~FacetProlongationTet() { ; }

    virtual void Update (const FESpace & fes){;} // no need to update 
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("FacetProlongationTet::CreateProlongationMatrix not implemented!");
    }
    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const;
    
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };
  
  // Facet Prolongation for FacetFESpace
  class FacetProlongationTet2 : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const meshTopology& et;
    public:
    FacetProlongationTet2(shared_ptr<MeshAccess> ama, const meshTopology & aet)
      : ma(ama), et(aet){;} 
    virtual ~FacetProlongationTet2() { ; }

    virtual void Update (const FESpace & fes){;} // no need to update 
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("FacetProlongationTet::CreateProlongationMatrix not implemented!");
    }
    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const;
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };
  

  /// Piecewise constant prolongaton.
  // Copied from ngsolve
  class ElementProlongation0 : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const meshTopology& et;
  public:
    ///
    ElementProlongation0(shared_ptr<MeshAccess> ama, const meshTopology & aet)
      : ma(ama), et(aet){;} 
    ///
    virtual ~ElementProlongation0() { ; }
  
    ///
    virtual void Update (const FESpace & fes) override
    { ; }

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override
    {
      FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

      int nc = et.GetCellLevel (finelevel-1);
      int nf = et.GetCellLevel (finelevel);

      for (int i = nc; i < nf; i++)
      {
        int parent = ma->GetParentElement (ElementId(VOL,i)).Nr();
        fv(i) = fv(parent);
      }

      for (int i = nf; i < fv.Size(); i++)
        fv(i) = 0;
    }

    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override
    {
      FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

      int nc = et.GetCellLevel (finelevel-1);
      int nf = et.GetCellLevel (finelevel);

      for (int i = nf-1; i >= nc; i--)
      {
        int parent = ma->GetParentElement (ElementId(VOL,i)).Nr();
        fv(parent) += fv(i);
        fv(i) = 0;
      }
    }
  };

}
