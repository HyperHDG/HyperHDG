#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_wave_ldgh.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include "parameters.hxx"
#include <map>

static const char help_msg[] = "convergence test regarding the wave equation\n";

PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, PetscReal* dat, PetscInt len) {
  PetscCall(PetscPrintf(com, msg));
  const PetscInt row_len = 10;
  for (PetscInt i = 0; i < len; i++) {
    if (i % row_len == 0)
       PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % .5e", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  return 0;
}

PetscErrorCode PetscPrin2i(MPI_Comm com, const char* msg, PetscInt* dat, PetscInt len) {
  PetscCall(PetscPrintf(com, msg));
  const PetscInt row_len = 10;
  for (PetscInt i = 0; i < len; i++) {
    if (i % row_len == 0)
      PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % 12d", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  return 0;
}

// must call VecRestoreSpan(x, span) after
PetscErrorCode VecGetSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p;
  PetscInt n;
  PetscCall(VecGetArray(x, &p));
  PetscCall(VecGetLocalSize(x, &n));
  span = {p, (size_t)n};
  return 0;
}

// must be called after each VecGetSpan(x, span)
PetscErrorCode VecRestoreSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p = span.data();
  PetscCall(VecRestoreArray(x, &p));
  span = std::span<PetscScalar>();
  return 0;
}

struct TestConfig {
  PetscReal tau;
  PetscReal theta;
  PetscInt iteration;
  PetscInt time_steps;
  PetscReal end_time;
  PetscReal dt;
};

template<int space_dim, int poly_deg>
PetscErrorCode test(TestConfig* cfg) {
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionWave<space_dim,poly_deg,2*poly_deg,TestWave,PetscReal>;
    using HDG = GlobalLoop::Hyperbolic<Top,Geo,NDes,LSol>;

    PetscInt N;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscFunctionBeginUser;
    HDG hdg((1 << cfg->iteration) * space_dim, {cfg->tau, cfg->theta, cfg->dt});

    zero_v = hdg.zero_vector();
    N = zero_v.size();
    temp = hdg.make_initial(zero_v);

    mat_coo = hdg.trace_to_flux_mat(0.);
    PetscCall(MatCreateFromOptions(PETSC_COMM_SELF, NULL, 1, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(MatSetPreallocationCOO(mat, mat_coo.value_vec.size(), (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, PETSC_TRUE));

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;
    PetscCall(VecGetSpan(rhs, rhs_span));
    PetscCall(VecGetSpan(sol, sol_span));

    for (PetscInt i = 0; i < cfg->time_steps; i++) {
        hdg.residual_flux2(std::span{zero_v}, rhs_span, (i+1)*cfg->dt);
        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, sol));
        hdg.set_data(sol_span, (i+1)*cfg->dt);
    }

    // zero_v unused
    temp2 = hdg.errors(zero_v, cfg->end_time);
    temp3 = hdg.norms(zero_v, cfg->end_time);
    for (size_t i = 0; i < temp3.size(); i++)
      temp3[i] = temp2[i] / temp3[i];

    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final abs error", temp2.data(), temp2.size()));
    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final rel error", temp3.data(), temp3.size()));

    PetscCall(VecRestoreSpan(rhs, rhs_span));
    PetscCall(VecRestoreSpan(sol, sol_span));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscFunctionReturn(0);
}

#define MAXLEN 100

int main(int argc, char **argv) {
    constexpr int space_dim = 1;
    constexpr int poly_deg = 3;
    PetscBool help = false, is_set;
    PetscReal err;
    TestConfig cfg;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscOptionsEnd();

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    cfg.tau = 1;
    cfg.theta = .25;
    cfg.iteration = 3;
    cfg.end_time = 1;
    cfg.time_steps = 100;
    cfg.dt = cfg.end_time / cfg.time_steps;

    test<space_dim,poly_deg>(&cfg);

    PetscCall(PetscFinalize());
    return 0;
}
