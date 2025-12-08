#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_wave_ldgh.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include "parameters.hxx"
#include <map>

static const char help_msg[] = "experiments regarding the wave equation\n";

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



int main(int argc, char **argv) {
    constexpr int space_dim = 1;
    constexpr int poly_deg = 3;
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionWave<space_dim,poly_deg,2*poly_deg,TestWave2,PetscReal>;
    using HDG = GlobalLoop::Hyperbolic<Top,Geo,NDes,LSol>;

    PetscBool help = false;
    PetscReal tau = 1; // HDG penalty
    PetscReal theta = .25; // one-step theta method
    PetscInt iteration = 1;
    PetscInt timesteps = 1;
    PetscReal end_time = 1;
    PetscReal dt = 0;

    PetscBool plot = true;
    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "wave";
    char plot_scale[PATH_MAX] = "0.95";

    PetscLogStage s_as, s_ts, s_rf;

    PetscBool is_set;
    PetscInt N;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_it = 0;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    std::vector<PetscInt> itemp;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsReal("-theta", "time-step averaging weight, 0 < theta <= 0.5, use theta=0.25 for CN", NULL, theta, &theta, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsInt("-i", "number of subdivisions of the domain, i >= 1", NULL, iteration, &iteration, &is_set));
    PetscCall(PetscOptionsInt("-ts", "number of timesteps", NULL, timesteps, &timesteps, &is_set));
    PetscCall(PetscOptionsReal("-T", "end time", NULL, end_time, &end_time, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-plot", "plot solution", NULL, plot, &plot, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    PetscCall(PetscLogStageRegister("Assembly", &s_as));
    PetscCall(PetscLogStageRegister("Timestepping", &s_ts));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));

    dt = end_time / timesteps;

    HDG hdg((1 << iteration) * space_dim, {tau, theta, dt});
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "true");
    hdg.plot_option("scale", plot_scale);

    zero_v = hdg.zero_vector();
    N = zero_v.size();
    temp = hdg.make_initial(zero_v);
    if (plot)
      hdg.plot_solution(temp, 0.);

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));

    PetscLogStagePush(s_as);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "assembly...\n"));
    mat_coo = hdg.trace_to_flux_mat(0.);
    PetscCall(MatCreateFromOptions(PETSC_COMM_SELF, NULL, 1, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(MatSetPreallocationCOO(mat, mat_coo.value_vec.size(), (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, PETSC_TRUE));
    PetscLogStagePop();

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(PetscPrintf(PETSC_COMM_SELF, "timestepping...\n"));
    PetscLogStagePush(s_ts);
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;
    PetscCall(VecGetSpan(rhs, rhs_span));
    PetscCall(VecGetSpan(sol, sol_span));
    for (PetscInt i = 0; i < timesteps; i++) {
        PetscLogStagePush(s_rf);
        hdg.residual_flux2(std::span{zero_v}, rhs_span, (i+1)*dt);

        PetscLogStagePop();
        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, sol));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        hdg.set_data(sol_span, (i+1)*dt);
        if (plot)
          hdg.plot_solution(sol_span, (i+1)*dt);
    }
    PetscLogStagePop();

    // zero_v unused
    temp2 = hdg.errors(zero_v, end_time);
    temp3 = hdg.norms(zero_v, end_time);
    for (size_t i = 0; i < temp3.size(); i++)
      temp3[i] = temp2[i] / temp3[i];
    avg_it = ((PetscReal)iterations) / timesteps;

    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final abs error", temp2.data(), temp2.size()));
    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final rel error", temp3.data(), temp3.size()));
    PetscCall(PetscPrin2i(PETSC_COMM_SELF, "tot iterations", &iterations, 1)); 
    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "avg iterations", &avg_it, 1)); 
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "wrote output to '%s/%s.*.vtu'\n", output_directory, output_filename));

    PetscCall(VecRestoreSpan(rhs, rhs_span));
    PetscCall(VecRestoreSpan(sol, sol_span));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}
