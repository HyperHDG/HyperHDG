// The module of tests/bindings.py, built by CMake (tests/CMakeLists.txt): LDG-H diffusion of
// degree 1 on the unit interval with GlobalLoop::Generic, bound as the class Diffusion.
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>

using Diffusion =
  GlobalLoop::Generic<Topology::Cubic<1, 1>,
                      Geometry::UnitCube<1, 1, double>,
                      NodeDescriptor::Cubic<1, 1>,
                      LocalSolver::Diffusion<1, 1, 2, tests::DiffusionSin, double> >;

NB_MODULE(hyperhdg_bindings_test, m)
{
  HyperHDG::bind_python<Diffusion, std::vector<unsigned int> >(m, "Diffusion");
}
