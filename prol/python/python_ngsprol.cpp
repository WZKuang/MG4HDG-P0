#include <core/utils.hpp>
#include <python_ngstd.hpp>

#include "../utils/python_utils.cpp"


PYBIND11_MODULE(ngsprol_py, m)
{
  cout << "importing ngs-prol" << NGSprol_VERSION << endl;
  ExportNgsx_utils(m);
}
