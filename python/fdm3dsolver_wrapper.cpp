#include <fdm3d/fdm3dsolver.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace {
namespace py = pybind11;

PYBIND11_MODULE(fdm3d, m) {
  py::class_<FDM3DSolver>(m, "FDM3D")
    .def(py::init<int, int, int>())
    .def("solve", &FDM3DSolver::solve)
    .def("getX", &FDM3DSolver::getX);
}
}  // namespace