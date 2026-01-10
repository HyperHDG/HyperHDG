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
};

template<int space_dim, int poly_deg>
PetscErrorCode test(TestConfig* cfg, PetscReal* err) {
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionWave<space_dim,poly_deg,2*poly_deg,TestWave3,PetscReal>;
    using HDG = GlobalLoop::Hyperbolic<Top,Geo,NDes,LSol>;

    PetscInt N;
    PetscReal dt = cfg->end_time / cfg->time_steps;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscFunctionBeginUser;
    HDG hdg((1 << cfg->iteration) * space_dim, {cfg->tau, cfg->theta, dt});

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
    PetscCall(KSPSetTolerances(ksp, 1e-16, 1e-16, PETSC_UNLIMITED, PETSC_UNLIMITED));
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;
    PetscCall(VecGetSpan(rhs, rhs_span));
    PetscCall(VecGetSpan(sol, sol_span));

    for (PetscInt i = 0; i < cfg->time_steps; i++) {
        hdg.residual_flux2(std::span{zero_v}, rhs_span, (i+1)*dt);
        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, sol));
        hdg.set_data(sol_span, (i+1)*dt);
    }

    // zero_v unused
    temp2 = hdg.errors(zero_v, cfg->end_time);
    temp3 = hdg.norms(zero_v, cfg->end_time);
    *err = temp2[0] / temp3[0];

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
    PetscReal errs[2] = {0}, alpha;
    TestConfig cfg;

    PetscInt n_it = MAXLEN;
    PetscInt n_ts = MAXLEN;
    PetscInt it[MAXLEN] = {0};
    PetscInt ts[MAXLEN] = {0};
    PetscReal theta = .25;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsIntArray("-i", "range of subdivisions of the domain, >= 1", NULL, it, &n_it, &is_set));
    PetscCall(PetscOptionsIntArray("-ts", "timestep exponents", NULL, ts, &n_ts, &is_set));
    PetscCall(PetscOptionsReal("-theta", "time-step averaging weight, 0 < theta <= 0.5, use theta=0.25 for CN", NULL, theta, &theta, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    cfg.tau = 1;
    cfg.theta = theta;
    cfg.end_time = 1;

    printf("# tau=%.5e\n", cfg.tau);
    printf("# theta=%.5e\n", cfg.theta);
    printf("# space_dim=%d\n", space_dim);
    printf("# poly_deg=%d\n", poly_deg);
    printf("i,ts,err,alpha\n");

    for (PetscInt i = 0; i < n_it; i++) {
      for (PetscInt j = 0; j < n_ts; j++) {
        cfg.iteration = it[i];
        cfg.time_steps = 1<<ts[j];

        std::swap(errs[0],errs[1]);
        test<space_dim,poly_deg>(&cfg, &errs[0]);
        alpha = std::log(errs[0]/errs[1]) / std::log(.5);
        printf("%04d,%04d,%.5e,%+.5e\n", 1<<cfg.iteration, cfg.time_steps, errs[0], alpha);
      }
    }

    PetscCall(PetscFinalize());
    return 0;
}
