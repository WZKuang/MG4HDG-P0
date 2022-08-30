#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>
#include <python_comp.hpp>
#include <bla.hpp>
#include <ngstd.hpp> // for Array
#include <nginterface.h>

using namespace ngsolve;
using namespace ngfem;

#include "../utils/NCprolongation.hpp"

using namespace ngcomp;

void ExportNgsx_utils(py::module &m)
{
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  
  
// meshTopology
  typedef shared_ptr<meshTopology> PymT;
  py::class_<meshTopology, PymT>
    (m, "meshTopology",
     docu_string(R"raw_string(mesh Topology for simplex mesh)raw_string"))
    .def("__init__",
        [](meshTopology *instance, shared_ptr<MeshAccess> ma, int dim){
        new (instance) meshTopology (ma, dim);
        },
        py::arg("mesh"), py::arg("dim")
        )
    .def("Update",
        [](shared_ptr<meshTopology> gt)
        {
        gt -> Update();
        }
        )
    .def("CreateSmoother",
        [](shared_ptr<meshTopology> gt, shared_ptr<BilinearForm> bfa, Flags &flags)
        {
        return gt -> CreateSmoother(bfa, flags);
        }
        )
    .def("CreateInnerSolve",
        [](shared_ptr<meshTopology> gt, shared_ptr<BilinearForm> bfa, Flags &flags, int finelevel)
        {
        return gt -> CreateInnerSolve(bfa, flags, finelevel);
        }
        );
  
  
  //Trig facet (FacetFESpace)
  typedef shared_ptr<FacetProlongationTrig> PyEdgePT3;
  py::class_<FacetProlongationTrig, PyEdgePT3, Prolongation>
    (m, "FacetProlongationTrig",
     docu_string(R"raw_string(
Facet Prolongation for Trig (P0)
)raw_string"))
    .def("__init__",
        [](FacetProlongationTrig *instance, shared_ptr<MeshAccess> ma,
          shared_ptr<meshTopology> mt)
        {
        new (instance) FacetProlongationTrig (ma, *mt);
        },
        py::arg("mesh"),
        py::arg("mesh Topology")
        )
    .def("GetInnerDofs",
        [](shared_ptr<FacetProlongationTrig> npt, int level)
        {
        auto inner = npt -> GetInnerDofs(level);
        return inner;
        },
        py::arg("level")
        )
    .def("Update",
        [](shared_ptr<FacetProlongationTrig> npt, shared_ptr<FESpace> fes)
        {
        npt -> Update(*fes);
        },
        py::arg("space")
        );
  

  //Trig facet (FacetFESpace)
  typedef shared_ptr<FacetProlongationTrig2> PyEdgePT32;
  py::class_<FacetProlongationTrig2, PyEdgePT32, Prolongation>
    (m, "FacetProlongationTrig2",
     docu_string(R"raw_string(
Facet Prolongation for Trig (NC-P0)
)raw_string"))
    .def("__init__",
        [](FacetProlongationTrig2 *instance, shared_ptr<MeshAccess> ma,
          shared_ptr<meshTopology> mt)
        {
        new (instance) FacetProlongationTrig2 (ma, *mt);
        },
        py::arg("mesh"),
        py::arg("mesh Topology")
        )
    .def("GetInnerDofs",
        [](shared_ptr<FacetProlongationTrig2> npt, int level)
        {
        auto inner = npt -> GetInnerDofs(level);
        return inner;
        },
        py::arg("level")
        )
    .def("Update",
        [](shared_ptr<FacetProlongationTrig2> npt, shared_ptr<FESpace> fes)
        {
        npt -> Update(*fes);
        },
        py::arg("space")
        );
  
  //Tet facet (FacetFESpace)
  typedef shared_ptr<FacetProlongationTet> PyFPT3;
  py::class_<FacetProlongationTet, PyFPT3, Prolongation>
    (m, "FacetProlongationTet",
     docu_string(R"raw_string(
Facet Prolongation for Tet (P0)
)raw_string"))
    .def("__init__",
        [](FacetProlongationTet *instance, shared_ptr<MeshAccess> ma,
          shared_ptr<meshTopology> mt)
        {
        new (instance) FacetProlongationTet (ma, *mt);
        },
        py::arg("mesh"),
        py::arg("mesh Topology")
        )
    .def("GetInnerDofs",
        [](shared_ptr<FacetProlongationTet> npt, int level)
        {
        auto inner = npt -> GetInnerDofs(level);
        return inner;
        },
        py::arg("level")
        )
    .def("Update",
        [](shared_ptr<FacetProlongationTet> npt, shared_ptr<FESpace> fes)
        {
        npt -> Update(*fes);
        },
        py::arg("space")
        );
  
  //Tet facet (nonconforming)
  typedef shared_ptr<FacetProlongationTet2> PyFPT32;
  py::class_<FacetProlongationTet2, PyFPT32, Prolongation>
    (m, "FacetProlongationTet2",
     docu_string(R"raw_string(
Facet Prolongation for Tet (NC)
)raw_string"))
    .def("__init__",
        [](FacetProlongationTet2 *instance, shared_ptr<MeshAccess> ma,
          shared_ptr<meshTopology> mt)
        {
        new (instance) FacetProlongationTet2 (ma, *mt);
        },
        py::arg("mesh"),
        py::arg("mesh Topology")
        )
    .def("GetInnerDofs",
        [](shared_ptr<FacetProlongationTet2> npt, int level)
        {
        auto inner = npt -> GetInnerDofs(level);
        return inner;
        },
        py::arg("level")
        )
    .def("Update",
        [](shared_ptr<FacetProlongationTet2> npt, shared_ptr<FESpace> fes)
        {
        npt -> Update(*fes);
        },
        py::arg("space")
        );
  

  //DG-P0
  typedef shared_ptr<ElementProlongation0> PyET;
  py::class_<ElementProlongation0, PyET, Prolongation>
    (m, "ElementProlongation0",
     docu_string(R"raw_string(Cell Prolongation)raw_string"))
    .def("__init__",
        [](ElementProlongation0 *instance, shared_ptr<MeshAccess> ma,
          shared_ptr<meshTopology> mt)
        {
        new (instance) ElementProlongation0 (ma, *mt);
        },
        py::arg("mesh"),
        py::arg("mesh Topology")
        );

}
