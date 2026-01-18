#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_wave1_ldgh.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>

#include "hdg_base.hxx"
#include "parameters.hxx"

static constexpr unsigned int poly_deg = 3;
template<unsigned int space_dim>
using HDGWave = GlobalLoop::Hyperbolic<
  Topology::Cubic<space_dim, space_dim>,
  Geometry::UnitCube<space_dim, space_dim, PetscReal>,
  NodeDescriptor::Cubic<space_dim, space_dim>,
  LocalSolver::DiffusionWave1<space_dim, poly_deg, 2*poly_deg, TestWave1, PetscReal>
>;
