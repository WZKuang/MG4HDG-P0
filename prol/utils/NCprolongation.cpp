#include "../utils/NCprolongation.hpp"
#include <algorithm>
//#include <bits/stdc++.h>
using namespace ngcomp;

// Function to sort character array b[]
// according to the order defined by a[]
void pairsort(int a[], float b[], int n)
{
    pair<int, float> pairt[n];
  
    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++) 
    {
        pairt[i].first = a[i];
        pairt[i].second = b[i];
    }
  
    // Sorting the pair array.
    sort(pairt, pairt + n);
      
    // Modifying original arrays
    for (int i = 0; i < n; i++) 
    {
        a[i] = pairt[i].first;
        b[i] = pairt[i].second;
    }
}

namespace ngmg
{
  // meshTopology (for trig)
  //
  void meshTopology :: Update()
  {
    size_t nv = ma->GetNV();
    size_t ned = ma->GetNEdges();
    size_t nfa = ma->GetNFaces();// # of faces (no intermediate faces)
    size_t ne = ma->GetNE();

    int level = ma->GetNLevels();
    if (level == nelevel.Size())
      return;
    size_t oldnv=0, oldned=0, oldnfa=0;
    Array<int> elnums;

    // # of vert/edge up to level 
    nvlevel.Append (nv); // all vertices <= level
    nelevel.Append (ned);// all edges <= level
    nflevel.Append (nfa);// all faces <= level
    nclevel.Append (ne);// all elements <= level

    // update edge nodes and node2edge hash-table 
    ClosedHashTable<INT<2>, int> node2edge(5*ned+10);
    edgenodes.SetSize0();
    for (size_t i = 0; i < ned; i++)
    {
      INT<2> edge = ma->GetEdgePNums (i);
      if (edge[0] > edge[1]) Swap (edge[0], edge[1]);
      node2edge.Set (edge, i); 
      edgenodes.Append (edge);
    }

    ClosedHashTable<INT<3>, int> node2face(5*nfa+10);
    facenodes.SetSize0();
    for (size_t i = 0; i < nfa; i++)
    {
      INT<3> face = ma->GetFacePNums (i);// always in accending order
      node2face.Set (face, i); 
      facenodes.Append (face);
    }

    shared_ptr<Array<int>> edge2dof = make_shared<Array<int>>(ned); 
    shared_ptr<Array<int>> face2dof = make_shared<Array<int>>(nfa); 
    // weights/coarse face info
    Array<float> cwfaces0;
    Array<int> ccfaces0;
    cwfaces0.SetSize(8);
    cwfaces0=0;
    ccfaces0.SetSize(8);
    ccfaces0=-1;

    // distribute vertex/edge levels
    if (level==1){// lowest level
      // calculate vertex levels
      vertexlevel.SetSize(nv);
      vertexlevel = 0; 

      if (dim==2){
        // calculate edge levels
        edgelevel.SetSize(ned);
        edgelevel = 0; 
        for (auto i : Range(ned))
          (*edge2dof)[i] = i;
        edgeDofNumber.Append(edge2dof);
        edgeIC.Append(0);// no inactive edges

        // update elemEdges
        Array<INT<3>> elem2edge, elem2vert;
        elem2edge.SetSize(ne);
        elem2vert.SetSize(ne);
        for (auto i: Range(ne)){
          auto edgenrs = ma->GetElEdges(i);
          auto vertnrs = ma->GetElVertices(i);
          elem2edge[i] = edgenrs;
          elem2vert[i] = vertnrs;
        }
        elemEdges.Append(elem2edge); // append
        elemVerts.Append(elem2vert); // append

        // update edgeElems
        Array<INT<2>> edge2elem, edge2lvl;
        edge2elem.SetSize(ned);
        edge2elem = -1;
        edge2lvl.SetSize(ned);
        edge2lvl = 0;
        for (auto i: Range(ned)){
          ma->GetEdgeElements (i, elnums);
          for (auto j=0; j< elnums.Size();j++)
            edge2elem[i][j] = elnums[j]; 
        }
        edgeElems.Append(edge2elem); // append
        edgeElemsLvl.Append(edge2lvl); // append

        ccedges.SetSize (ned); // coarsest edge no parents
        ccedges = -1;
      } else if (dim==3){
        // calculate edge levels
        facelevel.SetSize(nfa);
        facelevel = 0; 
        for (auto i : Range(nfa))
          (*face2dof)[i] = i;
        faceDofNumber.Append(face2dof);
        faceIC.Append(0);// no inactive edges

        // update elemEdges
        Array<INT<4>> elem2face, elem2vert;
        elem2face.SetSize(ne);
        elem2vert.SetSize(ne);
        for (auto i: Range(ne)){
          auto facenrs = ma->GetElFaces(i);
          auto vertnrs = ma->GetElVertices(i);
          elem2face[i] = facenrs;
          elem2vert[i] = vertnrs;
        }
        elemFaces.Append(elem2face); // append
        elemVerts3D.Append(elem2vert); // append

        // update edgeElems
        Array<INT<2>> facet2elem, facet2lvl;
        facet2elem.SetSize(nfa);
        facet2elem = -1;
        facet2lvl.SetSize(nfa);
        facet2lvl = 0;
        for (auto i: Range(nfa)){
          ma->GetFaceElements (i, elnums);
          for (auto j=0; j< elnums.Size();j++)
            facet2elem[i][j] = elnums[j]; 
        }
        faceElems.Append(facet2elem); // append
        faceElemsLvl.Append(facet2lvl); // append

        for (auto i:Range(nfa)){
          cwfaces.Append(cwfaces0);
          ccfaces.Append(ccfaces0);
        }
      }
    }
    else{
      // update vertex levels
      oldnv = vertexlevel.Size();
      vertexlevel0.SetSize(oldnv);
      vertexlevel0 = vertexlevel;
      vertexlevel.SetSize(nv);
      vertexlevel.Range (oldnv, nv) = level-1;
      vertexlevel.Range (0, oldnv) = vertexlevel0;

      // ************* BEGIN 2D edge parent nbrs
      // ************* BEGIN 2D edge parent nbrs
      // ************* BEGIN 2D edge parent nbrs
      if (dim==2){// update edge nbrs **necessary for hdg P0 prolongation
        // update edge levels
        oldned = edgelevel.Size();
        edgelevel0.SetSize(oldned);
        edgelevel0 = edgelevel;
        edgelevel.SetSize(ned);
        edgelevel.Range (oldned, ned) = level-1;
        edgelevel.Range (0, oldned) = edgelevel0;

        // update edge dofs
        auto iii = 0;
        auto jjj = ned-1;// OLD edge append to bottom 
        for (auto i : Range(ned))
        {
          ma->GetEdgeElements (i, elnums);
          if (elnums.Size()==0) {// this is an OLD edge, append to bottom
            (*edge2dof)[i] = jjj;
            jjj--;
          }
          else{// this is a NEW edge
            (*edge2dof)[i] = iii;
            iii++;
          }
        }
        edgeDofNumber.Append(edge2dof);
        edgeIC.Append(ned-iii);// inactive edges

        // update elemEdges
        Array<INT<3>> elem2edge, elem2vert;
        elem2edge.SetSize(ne);
        elem2vert.SetSize(ne);
        for (auto i: Range(ne)){
          auto edgenrs = ma->GetElEdges(i);
          auto vertnrs = ma->GetElVertices(i);
          elem2edge[i] = edgenrs;
          elem2vert[i] = vertnrs;
        }
        //cout << edgeelems<<endl;
        //cout << elem2edge<<endl;
        elemEdges.Append(elem2edge); // append
        elemVerts.Append(elem2vert); // append

        // update edgeElems
        //cout << edgeElems[level-2]<<endl;

        Array<INT<2>> edge2elem, edge2lvl;
        edge2elem.SetSize(ned);
        edge2elem=-1;
        edge2lvl.SetSize(ned);
        edge2lvl=-1;
        for (auto i: Range(ned)){
          ma->GetEdgeElements (i, elnums);
          if (elnums.Size()==0) {// old elem data
            edge2elem[i] = edgeElems[level-2][i];
            edge2lvl[i] = edgeElemsLvl[level-2][i];
          }else{
            for (auto j=0; j< elnums.Size();j++)
              edge2elem[i][j] = elnums[j]; 
            edge2lvl[i] = level-1;
          }
        }
        edgeElems.Append(edge2elem); // append
        edgeElemsLvl.Append(edge2lvl); // append
        //cout << edge2lvl << endl;
        //cout << edge2lvl<<endl;
        //cout << "*********"<<endl;

        // a new implementation
        // get coarse edge info...
        ccedges0.SetSize (oldned);
        ccedges0 = ccedges;
        ccedges.SetSize(ned);
        ccedges.Range(0, oldned) = ccedges0;
        ccedges.Range(oldned, ned)=-1;

        for (size_t i = oldned; i < ned; i++)
        {
          int pa[2], ii;
          int vv[4]={0,0,0,0};
          ii = 0;
          for (auto va: edgenodes[i]){
            if (vertexlevel[va] < level-1) {
              vv[ii++]=va;
              vv[ii++]=va;
            }else{
              ma->GetParentNodes (va, pa);
              vv[ii++] = pa[0]; 
              vv[ii++] = pa[1];
            }
          }
          //cout << vv[0]<<","<<vv[1]<<","<<vv[2]<<","<<vv[3]<<endl;
          // sort vv/get repeated
          sort(vv, vv+4);
          int ww[3]={-1,-1,-1}, wc[3]={0,0,0}, jj=0;
          ww[0]=vv[0]; wc[0]=1;
          for (int j=1;j<4;j++){
            if (vv[j]==ww[jj]) wc[jj]++;
            else {
              jj++;
              ww[jj] = vv[j];
              wc[jj]++;
            }
          }
          //for (int j=0;j<3;j++){
          //  cout <<ww[j]<<","<<wc[j]<<endl;
          //}
          INT<2> edge;
          if (ww[2] < 0){// bdry edge
            edge[0] = ww[0]; edge[1]=ww[1];
            auto paedgenr = node2edge.Get (edge);
            // FIXME lower lvl info
            // TODO: make an array of edgeelems
            auto elnums0 = edgeElems[level-2][paedgenr];
            auto lvl = edgeElemsLvl[level-2][paedgenr];
            int ii = 0;
            for (auto el: elnums0){
              if (el>-1){
                auto edgenrs = elemEdges[lvl[ii]][el];
                auto vertnrs = elemVerts[lvl[ii]][el];
                //cout << lvl[ii]<< ":"<<level-2<<endl;
                //auto edgenrs = elemEdges[level-2][el];
                //auto vertnrs = elemVerts[level-2][el];
                //cout << edgenrs<<":"<<vertnrs<<endl;
                for (auto j: Range(3)){
                  ccedges[i][j+6*ii] = edgenrs[j];
                  ccedges[i][j+3+6*ii] = 0;// set counter to zero
                }
                // determine locations of ww[0]/ww[1]
                int kk;
                for (auto j: Range(2)){
                  for (auto k:Range(3)){
                    if (vertnrs[k]==ww[j]){ kk=k; break;}
                  }
                  ccedges[i][3+kk+6*ii] = wc[j];
                }
                ii++;
              }
            }
          }else{// interior edge
            // FIXME THIS IS WRONG!!!! WHY
            //face[0] = ww[0]; face[1]=ww[1]; face[2]=ww[2];
            //auto el = node2face.Get (face);
            //cout << el<<endl; //FIXME
            edge[0] = ww[0]; edge[1]=ww[1];
            auto paedgenr = node2edge.Get (edge);
            auto elnums1 = edgeElems[level-2][paedgenr];
            auto lvl1 = edgeElemsLvl[level-2][paedgenr];
            int el,lvl;

            if (elnums1[1]==-1) {
              el=elnums1[0];
              lvl = lvl1[0];
            }
            else{
              edge[0] = ww[0]; edge[1]=ww[2];
              paedgenr = node2edge.Get (edge);
              auto elnums2 = edgeElems[level-2][paedgenr];
              auto lvl2 = edgeElemsLvl[level-2][paedgenr];
              if (elnums2[0]==elnums1[0] || elnums2[0]==elnums1[1]){
                el=elnums2[0]; lvl = lvl2[0];
              }
              else {
                el=elnums2[1]; lvl = lvl2[1];
              }
            }

            auto edgenrs = elemEdges[lvl][el];
            auto vertnrs = elemVerts[lvl][el];
            //cout << edgenrs<<":"<<vertnrs<<endl;
            for (auto j: Range(3)){
              ccedges[i][j] = edgenrs[j];
              ccedges[i][j+3] = 0;// set counter to zero
            }
            // determine locations of ww[0]/ww[1]
            int kk;
            for (auto j: Range(3)){
              for (auto k:Range(3)){
                if (vertnrs[k]==ww[j]){ kk=k; break;}
              }
              ccedges[i][3+kk] = wc[j];
            }
          }
        }
        //cout << ccedges<< endl;
      }
      else if (dim==3){// update face nbrs **necessary for hdg P0 prolongation
        // update edge levels
        oldnfa = facelevel.Size();
        facelevel0.SetSize(oldnfa);
        facelevel0 = facelevel;
        facelevel.SetSize(nfa);
        facelevel.Range (oldnfa, nfa) = level-1;
        facelevel.Range (0, oldnfa) = facelevel0;
        
        // update edge dofs
        auto iii = 0;
        auto jjj = nfa-1;// count backwards
        for (auto i : Range(nfa))
        {
          ma->GetFaceElements (i, elnums);
          if (elnums.Size()==0) {// this is an OLD face, append to bottom
            (*face2dof)[i] = jjj;
            jjj--;
          }
          else{// this is a NEW face
            (*face2dof)[i] = iii;
            iii++;
          }
        }
        faceDofNumber.Append(face2dof);
        faceIC.Append(nfa-iii);// inactive faces

        // update elemEdges
        Array<INT<4>> elem2facet, elem2vert;
        elem2facet.SetSize(ne);
        elem2vert.SetSize(ne);
        for (auto i: Range(ne)){
          auto facetnrs = ma->GetElFaces(i);
          auto vertnrs = ma->GetElVertices(i);
          elem2facet[i] = facetnrs;
          elem2vert[i] = vertnrs;
        }
        //cout << edgeelems<<endl;
        //cout << elem2edge<<endl;
        elemFaces.Append(elem2facet); // append
        elemVerts3D.Append(elem2vert); // append

        Array<INT<2>> facet2elem, facet2lvl;
        facet2elem.SetSize(nfa);
        facet2elem=-1;
        facet2lvl.SetSize(nfa);
        facet2lvl=-1;
        for (auto i: Range(nfa)){
          ma->GetFaceElements (i, elnums);
          if (elnums.Size()==0) {// old elem data
            facet2elem[i] = faceElems[level-2][i];
            facet2lvl[i] =  faceElemsLvl[level-2][i];
          }else{
            for (auto j=0; j< elnums.Size();j++)
              facet2elem[i][j] = elnums[j]; 
            facet2lvl[i] = level-1;
          }
        }
        faceElems.Append(facet2elem); // append
        faceElemsLvl.Append(facet2lvl); // append

        // a new implementation
        // get coarse face info...
        for (size_t i = oldnfa; i < nfa; i++)
        {
          cwfaces0 = 0;// initialize weights
          ccfaces0 = -1;// initialize weights
          int pa[2], pa2[2], ii, jj;
          int vv[8]={0,0,0,0,0,0,-1,-1};
          float vw[8]={1,1,1,1,1,1,0,0};
          ii = 0; jj =6;
          for (auto va: facenodes[i]){
            if (vertexlevel[va] < level-1) {
              vv[ii++]=va;
              vv[ii++]=va;
            }else{
              ma->GetParentNodes (va, pa);
              vv[ii++] = pa[0]; 
              vv[ii++] = pa[1];
              // do something special if parent node is NEW
              if (vertexlevel[pa[1]] == level-1){
                ma->GetParentNodes (pa[1], pa2);
                vv[ii-1] = pa2[0];
                vw[ii-1] = 0.5;
                vv[jj]=pa2[1];
                vw[jj] = 0.5;
                jj++;
              }
            }
          }
          pairsort(vv, vw, 8);
          // sort vv/get repeated
          int ww[4]={-1,-1,-1,-1};
          float wc[4]={0,0,0,0};
          jj = 0;
          int flag = 0;
          // TODO
          for (int j=0;j<8;j++){
            if (vv[j]<0) continue;
            // initial value
            if (ww[jj]==-1 && flag==0) {
              ww[jj] = vv[j];
              flag = 1;
            }
            if (vv[j]==ww[jj]) wc[jj] += vw[j];
            else {
              jj++;
              ww[jj] = vv[j];
              wc[jj] += vw[j];
            }
          }
          INT<3> face;
          if (ww[3] < 0){// bdry face
            face[0] = ww[0]; face[1]=ww[1]; face[2] =ww[2];
            auto pafacenr = node2face.Get (face);
            auto elnums0 = faceElems[level-2][pafacenr];
            auto lvl = faceElemsLvl[level-2][pafacenr];
            int ii = 0;
            for (auto el: elnums0){
              if (el>-1){
                auto facenrs = elemFaces[lvl[ii]][el];
                auto vertnrs = elemVerts3D[lvl[ii]][el];
                for (auto j: Range(4)){
                  //ccfaces[i][j+4*ii] = facenrs[j];
                  ccfaces0[j+4*ii] = facenrs[j];
                }
                // determine locations of ww[0]/ww[1]/ww[2]
                int kk;
                for (auto j: Range(3)){
                  for (auto k:Range(4)){
                    if (vertnrs[k]==ww[j]){ kk=k; break;}
                  }
                  cwfaces0[kk+4*ii] = wc[j]; // weights 
                }
                ii++;
              }
            }
          }
          else{// interior face
            face[0] = ww[0]; face[1]=ww[1]; face[2] = ww[2];
            auto pafacenr = node2face.Get (face);
            auto elnums1 = faceElems[level-2][pafacenr];
            auto lvl1 = faceElemsLvl[level-2][pafacenr];
            int el,lvl;

            if (elnums1[1]==-1) {
              el=elnums1[0];
              lvl = lvl1[0];
            }
            else{
              face[2] = ww[3];
              pafacenr = node2face.Get (face);
              auto elnums2 = faceElems[level-2][pafacenr];
              auto lvl2 = faceElemsLvl[level-2][pafacenr];
              if (elnums2[0]==elnums1[0] || elnums2[0]==elnums1[1]){
                el=elnums2[0]; lvl = lvl2[0];
              }
              else {
                el=elnums2[1]; lvl = lvl2[1];
              }
            }

            auto facenrs = elemFaces[lvl][el];
            auto vertnrs = elemVerts3D[lvl][el];
            for (auto j: Range(4)){
              //ccfaces[i][j] = facenrs[j];
              ccfaces0[j] = facenrs[j];
            }
            // determine locations of ww[0]/ww[1]
            int kk;
            for (auto j: Range(4)){
              for (auto k:Range(4)){
                if (vertnrs[k]==ww[j]){ kk=k; break;}
              }
              cwfaces0[kk] = wc[j]; // weights 
            }
          }
          cwfaces.Append(cwfaces0);
          ccfaces.Append(ccfaces0);
        }
      }
    }
    // check labelling error
    //(*testout) << "vert level: " << vertexlevel << endl;
  }
  
  shared_ptr<BaseBlockJacobiPrecond> meshTopology ::CreateSmoother (shared_ptr<BilinearForm> bfa, 
        Flags & flags) const
  {
    // new: blocktypes, specified in fespace
    if (bfa->UsesEliminateInternal()){
      flags.SetFlag("eliminate_internal");
    }

    //shared_ptr<Table<int>> blocks = bfa->GetFESpace()->CreateSmoothingBlocks(flags);
    // NO PARALLEL FIXME
    //return dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
    //  .CreateBlockJacobiPrecond(blocks, 0, true, bfa->GetFESpace()->GetFreeDofs());
    //;
    
    // FIXME: blocking 
    auto fes = bfa->GetFESpace();
    auto freedofs = fes->GetFreeDofs(true);

    FilteredTableCreator creator(freedofs.get());
    
    if (flags.GetStringFlag("blocktype") == "edgepatch")
      {
        
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            // EFI
            if (dim==2)
              for (size_t i : Range(ma->GetNEdges()))        
              {
                // Ng_Node<1> edge = ma->GetNode<1> (i);

                fes->GetDofNrs (NodeId(NT_EDGE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    creator.Add (i, d);
              }
            
            if (dim==3)
              for (size_t i : Range(ma->GetNFaces()))        
              {
                // Ng_Node<2> face = ma->GetNode<2> (i);

                fes->GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    for (auto e : ma->GetFaceEdges(i))
                      creator.Add (e, d);
              }
            
            //if(ma->GetDimension() == 3)
            //  for(size_t i : Range(ma->GetNE()))
            //    {
            //      fes->GetDofNrs(NodeId(NT_CELL, i), dofs);
            //      auto eledges = ma->GetElEdges(ElementId(VOL, i));
            //      for(auto d : dofs)
            //        if(fes->IsRegularDof(d))
            //          for(auto v : eledges)
            //            creator.Add(v, d);
            //    }
          }
      }
    else if (flags.GetStringFlag("blocktype") == "point")
      {
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            if (dim==2) // 2D
              for (size_t i : Range(ma->GetNEdges()))
              {
                Ng_Node<1> edge = ma->GetNode<1> (i);

                fes->GetDofNrs (NodeId(NT_EDGE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    creator.Add (i, d);
              }

            if (dim==3) 
              for (size_t i : Range(ma->GetNFaces()))        
              {
                Ng_Node<2> face = ma->GetNode<2> (i);

                fes->GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    creator.Add (i, d);
              }
          }
      }
    else
      { // default is vertexpatch
        Array<DofId> dofs;
        for ( ; !creator.Done(); creator++)
          {
            // VEFI
            //for (size_t i : Range(ma->GetNV()))
            //  {
            //    fes->GetDofNrs (NodeId(NT_VERTEX, i), dofs);
            //    for (auto d : dofs)
            //      if (d>=0)
            //        creator.Add (i, d);
            //  }
            if (dim==2) // 2D
              for (size_t i : Range(ma->GetNEdges()))        
              {
                Ng_Node<1> edge = ma->GetNode<1> (i);

                fes->GetDofNrs (NodeId(NT_EDGE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    for (int k = 0; k < 2; k++)
                      creator.Add (edge.vertices[k], d);
              }
            if (dim==3) 
              for (size_t i : Range(ma->GetNFaces()))        
              {
                Ng_Node<2> face = ma->GetNode<2> (i);

                fes->GetDofNrs (NodeId(NT_FACE, i), dofs);
                for (auto d : dofs)
                  if (d>=0)
                    for (int k = 0; k < face.vertices.Size(); k++)
                      creator.Add (face.vertices[k], d);
              }
            
            //if(ma->GetDimension() == 3)
            //  for(size_t i : Range(ma->GetNE()))
            //    {
            //      fes->GetDofNrs(NodeId(NT_CELL, i), dofs);
            //      auto elverts = ma->GetElVertices(ElementId(VOL, i));
            //      for(auto d : dofs)
            //      if (fes->IsRegularDof(d))
            //          for(auto v : elverts)
            //            creator.Add(v, d);
            //    }
          }
      }

    // return make_shared<Table<int>> (creator.MoveTable());
    Table<int> table = creator.MoveTable();
    // NO PARALLEL WHAT IS PARALLEL??? FIXME
    return dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
      .CreateBlockJacobiPrecond(make_shared<Table<int>> (move(table)),
          0, true, bfa->GetFESpace()->GetFreeDofs());
    ;
  }
  
  
  // TODO: get inner dofs
  shared_ptr<BaseBlockJacobiPrecond> meshTopology ::CreateInnerSolve (shared_ptr<BilinearForm> bfa, 
        Flags &flags, int finelevel) const
  {
    throw Exception("meshTopology::InnerSolve not implemented!");
  }

  // 2D trig
  // Facet (FacetFE)
  shared_ptr<BitArray> FacetProlongationTrig ::  GetInnerDofs (int finelevel) const
  {
    size_t nc = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t nf = et.GetEdgeLevel(finelevel);
    auto edgeDofsF = et.GetEdgeDofNumber(finelevel);// fine edge dofs
    size_t nf0 = et.GetEdgeIC(finelevel);
    auto nfB = nf-nf0;
    //cout << nc <<","<<nf<<endl;

    BitArray inner(nfB);
    inner.Clear();
    for (size_t i = nc; i < nf; i++)
    {
      auto ccedges = et.CCEdges(i);
      if (ccedges[3] != 0 and ccedges[4] != 0 and ccedges[5] != 0)// this is an interior edge
            inner.SetBit((*edgeDofsF)[i]);
    }
    //cout << inner << endl;
    return make_shared<BitArray> (inner);
  }
  
  void FacetProlongationTrig ::  ProlongateInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t neF = et.GetEdgeLevel(finelevel);
    auto edgeDofsC = et.GetEdgeDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetEdgeDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetEdgeIC(finelevel-1);
    size_t neF0 = et.GetEdgeIC(finelevel);
    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    // temporary vector
    shared_ptr<BaseVector> tmp_vec = make_shared<VVector<double>>(neF);// all edges
    FlatVector<> fw = tmp_vec->FV<double>();

    // STEP 1: fv-->fw (copy active edge data)
    fw = 0;
    auto jjj = 0;
    for (auto i = 0; i < neC; i++){
      if ((*edgeDofsC)[i]<neCB)
        fw((*edgeDofsF)[i]) = fv(jjj++);
    }
    
    // STEP 2: fw-->fw update
    for (auto i = neC; i < neF; i++)
    {
      
      // NEW implementation
      auto ccedges = et.CCEdges(i);
      double fac = 0.25, e0, e1, e2, vv[3];
      int els = 1;
      if (ccedges[6]>-1) {
        fac = 0.125; els=2;
      }
      
      for (auto i0: Range(els)){
        e0 = fw((*edgeDofsF)[ccedges[0+6*i0]]);
        e1 = fw((*edgeDofsF)[ccedges[1+6*i0]]);
        e2 = fw((*edgeDofsF)[ccedges[2+6*i0]]);
        vv[0] = e2+e0-e1;
        vv[1] = e2+e1-e0;
        vv[2] = e1+e0-e2;
        for (auto j: Range(3)){
          fw((*edgeDofsF)[i]) += fac*vv[j]*ccedges[3+j+6*i0];
        }
      }
    }

    // STEP 3: fw-->fv update
    fv=0;
    jjj = 0;
    for (auto i = 0; i < neF; i++)
    {
      if ((*edgeDofsF)[i]<neFB)
        fv(jjj++) = fw((*edgeDofsF)[i]);
    }
  }

  void FacetProlongationTrig::RestrictInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t neF = et.GetEdgeLevel(finelevel);
    auto edgeDofsC = et.GetEdgeDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetEdgeDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetEdgeIC(finelevel-1);
    size_t neF0 = et.GetEdgeIC(finelevel);
    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    // temporary vector
    shared_ptr<BaseVector> tmp_vec = make_shared<VVector<double>>(neF);// all edges
    FlatVector<> fw = tmp_vec->FV<double>();

    // STEP 3: fw-->fv update
    fw=0;
    auto jjj = 0;
    for (auto i = 0; i < neF; i++)
    {
      if ((*edgeDofsF)[i]<neFB)
        fw((*edgeDofsF)[i]) = fv(jjj++);
    }
    fv=0;

    // STEP 2: fw-->fw update
    for (auto i = neF-1; i >= neC; i--)
    {
      // NEW IMPLEMENTATION
      auto ccedges = et.CCEdges(i);
      double fac = 1, e0, e1, e2, vv[3], val;
      int els = 1;
      if (ccedges[6]>-1) {
        fac = 0.5; els=2;
      }

      for (auto i0: Range(els)){
        vv[0]=0; vv[1]=0; vv[2]=0;
        // update vert values
        for (auto j: Range(3)){
          // pay attention to weights
          vv[j] = fac*fw((*edgeDofsF)[i])*(ccedges[3+j+6*i0]-1);
        }
        // v--> e
        fw((*edgeDofsF)[ccedges[0+6*i0]]) += 0.5*(vv[0]+vv[2]);
        fw((*edgeDofsF)[ccedges[1+6*i0]]) += 0.5*(vv[1]+vv[2]);
        fw((*edgeDofsF)[ccedges[2+6*i0]]) += 0.5*(vv[0]+vv[1]);
      }
    }

    // STEP 1: fv-->fw (copy active edge data)
    jjj = 0;
    for (auto i = 0; i < neC; i++){
      if ((*edgeDofsC)[i]<neCB)
        fv(jjj++) = fw((*edgeDofsF)[i]);
    }
  }


  // 2D trig,NC
  // Facet (FacetFE)
  shared_ptr<BitArray> FacetProlongationTrig2 ::  GetInnerDofs (int finelevel) const
  {
    size_t nc = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t nf = et.GetEdgeLevel(finelevel);
    
    BitArray inner(nf);
    inner.Clear();
    for (size_t i = nc; i < nf; i++)
    {
      auto ccedges = et.CCEdges(i);
      if (ccedges[3] != 0 and ccedges[4] != 0 and ccedges[5] != 0)// this is an interior edge
            inner.SetBit(i);
    }
    return make_shared<BitArray> (inner);
  }

  void FacetProlongationTrig2 ::  ProlongateInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t neF = et.GetEdgeLevel(finelevel);
    auto edgeDofsC = et.GetEdgeDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetEdgeDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetEdgeIC(finelevel-1);
    size_t neF0 = et.GetEdgeIC(finelevel);
    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    fv.Range(neC, fv.Size()) = 0;
    
    // STEP 2: fw-->fw update
    for (auto i = neC; i < neF; i++)
    {
      
      // NEW implementation
      auto ccedges = et.CCEdges(i);
      double fac = 0.25, e0, e1, e2, vv[3];
      int els = 1;
      if (ccedges[6]>-1) {
        fac = 0.125; els=2;
      }
      
      for (auto i0: Range(els)){
        e0 = fv(ccedges[0+6*i0]);
        e1 = fv(ccedges[1+6*i0]);
        e2 = fv(ccedges[2+6*i0]);
        vv[0] = e2+e0-e1;
        vv[1] = e2+e1-e0;
        vv[2] = e1+e0-e2;
        for (auto j: Range(3)){
          fv(i) += fac*vv[j]*ccedges[3+j+6*i0];
        }
      }
    }
    
    // FIXME: remove data on inactive edges
    for (size_t i = 0; i < neF; i++)
      if ((*edgeDofsF)[i]>=neFB)
        fv(i) = 0;
  }

  void FacetProlongationTrig2::RestrictInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetEdgeLevel(finelevel-1);//-1 level
    size_t neF = et.GetEdgeLevel(finelevel);
    auto edgeDofsC = et.GetEdgeDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetEdgeDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetEdgeIC(finelevel-1);
    size_t neF0 = et.GetEdgeIC(finelevel);
    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    fv.Range(neF, fv.Size()) = 0;

    // remove coarse grid data
    for (size_t i = 0; i < neF; i++)
      if ((*edgeDofsF)[i]>=neFB)
        fv(i) = 0;

    // STEP 2: fw-->fw update
    for (auto i = neF-1; i >= neC; i--)
    {
      // NEW IMPLEMENTATION
      auto ccedges = et.CCEdges(i);
      double fac = 1, e0, e1, e2, vv[3], val;
      int els = 1;
      if (ccedges[6]>-1) {
        fac = 0.5; els=2;
      }

      for (auto i0: Range(els)){
        vv[0]=0; vv[1]=0; vv[2]=0;
        // update vert values
        for (auto j: Range(3)){
          // pay attention to weights
          vv[j] = fac*fv(i)*(ccedges[3+j+6*i0]-1);
        }
        // v--> e
        fv(ccedges[0+6*i0]) += 0.5*(vv[0]+vv[2]);
        fv(ccedges[1+6*i0]) += 0.5*(vv[1]+vv[2]);
        fv(ccedges[2+6*i0]) += 0.5*(vv[0]+vv[1]);
      }
    }
  }







  // 3D tet
  shared_ptr<BitArray> FacetProlongationTet ::  GetInnerDofs (int finelevel) const
  {
    size_t nc = et.GetFaceLevel(finelevel-1);//-1 level
    size_t nf = et.GetFaceLevel(finelevel);
    auto edgeDofsF = et.GetFaceDofNumber(finelevel);// fine edge dofs
    size_t nf0 = et.GetFaceIC(finelevel);
    auto nfB = nf-nf0;

    BitArray inner(nfB);
    inner.Clear();
    for (size_t i = nc; i < nf; i++)
    {
      if (et.CWFaces(i,0) !=0 && et.CWFaces(i,1) !=0 && 
          et.CWFaces(i,2) !=0 && et.CWFaces(i,3) !=0)
            inner.SetBit((*edgeDofsF)[i]);
    }
    return make_shared<BitArray> (inner);
  }
  
  void FacetProlongationTet ::  ProlongateInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetFaceLevel(finelevel-1);//-1 level
    size_t neF = et.GetFaceLevel(finelevel);
    auto edgeDofsC = et.GetFaceDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetFaceDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetFaceIC(finelevel-1); // in active faces
    size_t neF0 = et.GetFaceIC(finelevel);

    auto neCB = neC-neC0;
    auto neFB = neF-neF0;
    
    FlatVector<> fv = v.FV<double>();
    // temporary vector
    shared_ptr<BaseVector> tmp_vec = make_shared<VVector<double>>(neF);// all edges
    FlatVector<> fw = tmp_vec->FV<double>();
    // STEP 1: fv-->fw (copy active edge data)
    fw = 0;
    auto jjj = 0;
    for (auto i = 0; i < neC; i++){
      if ((*edgeDofsC)[i]<neCB)
        fw((*edgeDofsF)[i]) = fv(jjj++);
    }
    
    // STEP 2: fw-->fw update
    for (auto i = neC; i < neF; i++)
    {
      
      // NEW implementation
      double fac = 1.0/6.0, e0, e1, e2, e3, ee, vv[4];
      int els = 1;
      if (et.CCFaces(i, 4)>-1) {
        fac /= 2; els=2;
      }
      
      for (auto i0: Range(els)){
        e0 = fw((*edgeDofsF)[et.CCFaces(i, 0+4*i0)]);
        e1 = fw((*edgeDofsF)[et.CCFaces(i, 1+4*i0)]);
        e2 = fw((*edgeDofsF)[et.CCFaces(i, 2+4*i0)]);
        e3 = fw((*edgeDofsF)[et.CCFaces(i, 3+4*i0)]);
        ee = e0+e1+e2+e3;
        vv[0] = ee-3*e0;
        vv[1] = ee-3*e1;
        vv[2] = ee-3*e2;
        vv[3] = ee-3*e3;
        for (auto j: Range(4)){
          fw((*edgeDofsF)[i]) += fac*vv[j]*et.CWFaces(i,j+4*i0);
        }
      }
    }

    // STEP 3: fw-->fv update
    fv=0;
    jjj = 0;
    for (auto i = 0; i < neF; i++)
    {
      if ((*edgeDofsF)[i]<neFB)
        fv(jjj++) = fw((*edgeDofsF)[i]);
    }
  }

  void FacetProlongationTet::RestrictInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetFaceLevel(finelevel-1);//-1 level
    size_t neF = et.GetFaceLevel(finelevel);
    auto edgeDofsC = et.GetFaceDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetFaceDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetFaceIC(finelevel-1); // in active faces
    size_t neF0 = et.GetFaceIC(finelevel);

    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    // temporary vector
    shared_ptr<BaseVector> tmp_vec = make_shared<VVector<double>>(neF);// all edges
    FlatVector<> fw = tmp_vec->FV<double>();

    // STEP 3: fw-->fv update
    fw=0;
    auto jjj = 0;
    for (auto i = 0; i < neF; i++)
    {
      if ((*edgeDofsF)[i]<neFB)
        fw((*edgeDofsF)[i]) = fv(jjj++);
    }
    fv=0;
    
    // STEP 2: fw-->fw update
    for (auto i = neF-1; i >= neC; i--)
    {
      // NEW IMPLEMENTATION
      double fac0 = 1, e0, e1, e2, e3, vv[4];
      int els = 1;
      if (et.CCFaces(i, 4)>-1) {
        fac0 = 0.5; els=2;
      }
      //cout << ccedges << endl;
      for (auto i0: Range(els)){
        vv[0]=0; vv[1]=0; vv[2]=0; vv[3] = 0;
        // update vert values
        for (auto j: Range(4)){
          // pay attention to weights
          vv[j] = fac0*fw((*edgeDofsF)[i])*(et.CWFaces(i,j+4*i0)*1.5-2);
        }
        fw((*edgeDofsF)[et.CCFaces(i, 0+4*i0)]) += (vv[1]+vv[2]+vv[3])/3;
        fw((*edgeDofsF)[et.CCFaces(i, 1+4*i0)]) += (vv[0]+vv[2]+vv[3])/3;
        fw((*edgeDofsF)[et.CCFaces(i, 2+4*i0)]) += (vv[0]+vv[1]+vv[3])/3;
        fw((*edgeDofsF)[et.CCFaces(i, 3+4*i0)]) += (vv[0]+vv[1]+vv[2])/3;
      }
    }

    // STEP 1: fv-->fw (copy active edge data)
    jjj = 0;
    for (auto i = 0; i < neC; i++){
      if ((*edgeDofsC)[i]<neCB)
        fv(jjj++) = fw((*edgeDofsF)[i]);
    }
  }
  

  // 3D tet//NC
  shared_ptr<BitArray> FacetProlongationTet2 ::  GetInnerDofs (int finelevel) const
  {
    size_t nc = et.GetFaceLevel(finelevel-1);//-1 level
    size_t nf = et.GetFaceLevel(finelevel);

    BitArray inner(nf);
    inner.Clear();
    for (size_t i = nc; i < nf; i++)
    {
      if (et.CWFaces(i,0) !=0 && et.CWFaces(i,1) !=0 && 
          et.CWFaces(i,2) !=0 && et.CWFaces(i,3) !=0)
            inner.SetBit(i);
    }
    return make_shared<BitArray> (inner);
  }
  
  void FacetProlongationTet2 ::  ProlongateInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetFaceLevel(finelevel-1);//-1 level
    size_t neF = et.GetFaceLevel(finelevel);
    auto edgeDofsC = et.GetFaceDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetFaceDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetFaceIC(finelevel-1); // in active faces
    size_t neF0 = et.GetFaceIC(finelevel);

    auto neCB = neC-neC0;
    auto neFB = neF-neF0;
    
    FlatVector<> fv = v.FV<double>();
    fv.Range(neC, fv.Size()) = 0;
    
    // STEP 2: fw-->fw update
    for (auto i = neC; i < neF; i++)
    {
      
      // NEW implementation
      double fac = 1.0/6.0, e0, e1, e2, e3, ee, vv[4];
      int els = 1;
      if (et.CCFaces(i, 4)>-1) {
        fac /= 2; els=2;
      }
      
      for (auto i0: Range(els)){
        e0 = fv(et.CCFaces(i, 0+4*i0));
        e1 = fv(et.CCFaces(i, 1+4*i0));
        e2 = fv(et.CCFaces(i, 2+4*i0));
        e3 = fv(et.CCFaces(i, 3+4*i0));
        ee = e0+e1+e2+e3;
        vv[0] = ee-3*e0;
        vv[1] = ee-3*e1;
        vv[2] = ee-3*e2;
        vv[3] = ee-3*e3;
        for (auto j: Range(4)){
          fv(i) += fac*vv[j]*et.CWFaces(i,j+4*i0);
        }
      }
    }
    
    // FIXME: remove data on inactive edges
    for (size_t i = 0; i < neF; i++)
      if ((*edgeDofsF)[i]>=neFB)
        fv(i) = 0;

  }

  void FacetProlongationTet2::RestrictInline (int finelevel, BaseVector & v) const
  {
    size_t neC = et.GetFaceLevel(finelevel-1);//-1 level
    size_t neF = et.GetFaceLevel(finelevel);
    auto edgeDofsC = et.GetFaceDofNumber(finelevel-1);// coarse edge dofs
    auto edgeDofsF = et.GetFaceDofNumber(finelevel);// fine edge dofs
    size_t neC0 = et.GetFaceIC(finelevel-1); // in active faces
    size_t neF0 = et.GetFaceIC(finelevel);

    auto neCB = neC-neC0;
    auto neFB = neF-neF0;

    FlatVector<> fv = v.FV<double>();
    fv.Range(neF, fv.Size()) = 0;

    // remove coarse grid data
    for (size_t i = 0; i < neF; i++)
      if ((*edgeDofsF)[i]>=neFB)
        fv(i) = 0;
    
    // STEP 2: fw-->fw update
    for (auto i = neF-1; i >= neC; i--)
    {
      // NEW IMPLEMENTATION
      double fac0 = 1, e0, e1, e2, e3, vv[4];
      int els = 1;
      if (et.CCFaces(i, 4)>-1) {
        fac0 = 0.5; els=2;
      }
      //cout << ccedges << endl;
      for (auto i0: Range(els)){
        vv[0]=0; vv[1]=0; vv[2]=0; vv[3] = 0;
        // update vert values
        for (auto j: Range(4)){
          // pay attention to weights
          vv[j] = fac0*fv(i)*(et.CWFaces(i,j+4*i0)*1.5-2);
        }
        fv(et.CCFaces(i, 0+4*i0)) += (vv[1]+vv[2]+vv[3])/3;
        fv(et.CCFaces(i, 1+4*i0)) += (vv[0]+vv[2]+vv[3])/3;
        fv(et.CCFaces(i, 2+4*i0)) += (vv[0]+vv[1]+vv[3])/3;
        fv(et.CCFaces(i, 3+4*i0)) += (vv[0]+vv[1]+vv[2])/3;
      }
    }
  }


}

