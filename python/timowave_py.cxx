// Explicitly instantiated nanobind module for the timowave/wave4 path: this file spells out
// the (poly_deg, stages) combinations that get compiled -- add a bind_python line to compile
// another one. Serial, PETSc-free (real stage-1 systems are plain double, multi-stage Gauss
// representatives are std::complex<double>); python consumes the COO system with scipy.
//
// Mirrors the HDGTimoWave alias of experiments/timowave.cxx with double instead of PetscReal.

// first include: Python.h (via nanobind) must precede any libc feature-test macros
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <HyperHDG/global_loop/hyperbolic.hxx>
#include <HyperHDG/local_solver/timowave.hxx>

#include "timowave4.hxx"

#include <complex>

template <unsigned int poly_deg, unsigned int stages, typename ScalarT>
using HDGTimoWave = GlobalLoop::Hyperbolic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2 * poly_deg, TestTimoWave4, double, stages>,
  std::vector<ScalarT>>;

NB_MODULE(timowave_py, m)
{
  using Complex = std::complex<double>;

  // wave4 spatial phase, the -wave4_px option of the PETSc driver (TestTimoWave4::Init);
  // set it before constructing a solver.
  m.def(
    "wave4_set_px", [](double px) { TestTimoWave4<3, double>::px = px; },
    "spatial phase px of the wave4 standing-wave factor cos(w*(x+y+z) + px)");

  HyperHDG::bind_python<HDGTimoWave<1, 1, double>>(m, "TimoWave4_P1S1");
  HyperHDG::bind_python<HDGTimoWave<2, 2, Complex>>(m, "TimoWave4_P2S2");
}
