// The nanobind module of examples/diffusion_elliptic_nb_precompiled.py, built by CMake
// (examples/CMakeLists.txt) as diffusion_elliptic_nb: LDG-H diffusion of degree 1 on the unit
// square with GlobalLoop::Generic, bound as the class Diffusion.
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <examples/parameters/diffusion.hxx>

using Diffusion =
  GlobalLoop::Generic<Topology::Cubic<2, 2>,
                      Geometry::UnitCube<2, 2, double>,
                      NodeDescriptor::Cubic<2, 2>,
                      LocalSolver::Diffusion<2, 1, 2, HG<2>::DiffusionElliptic, double> >;

NB_MODULE(diffusion_elliptic_nb, m)
{
  HyperHDG::bind_python<Diffusion, std::vector<unsigned int> >(m, "Diffusion");
}
