#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>

#include <HyperHDG/local_solver/timowave.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>

#include "parameters.hxx"
#include "hdg_base.hxx"
#include "prin2.hxx"
#include "net2as.hxx"

static const char help_msg[] = "experiments regarding the wave equation\n";

template<unsigned int poly_deg, template<unsigned int, typename param_float_t> typename Test>
using HDGTimoWave = GlobalLoop::Hyperbolic<
  Topology::File<1,3>,
  Geometry::File<1,3>,
  NodeDescriptor::File<1,3>,
  LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2*poly_deg, Test, PetscReal>
>;

// hdg must be deallocated with `delete`
PetscErrorCode PetscHDGCreate(
    PetscInt poly_deg, PetscInt test,
    const char *path, PetscReal tau, PetscReal theta, PetscReal dt,
    HDGBase **hdg
) {
  int i = poly_deg*10 + test;
  switch(i) {
  case 14: *hdg = new HDGWrapper(HDGTimoWave<1,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 24: *hdg = new HDGWrapper(HDGTimoWave<2,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 34: *hdg = new HDGWrapper(HDGTimoWave<3,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 30: *hdg = new HDGWrapper(HDGTimoWave<3,TimoWaveClamped>(path, {tau, theta, dt})); return 0;
  case 64: *hdg = new HDGWrapper(HDGTimoWave<6,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  default:
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "unsupported: poly_deg = %d, test = %d", poly_deg, test);
  }

  return 0;
}

int main(int argc, char **argv) {
    PetscBool help = false, is_set;
    PetscInt nt = 1, nx = 1, poly_deg = 3;
    PetscInt N;            // global system size
    PetscReal tau = 1;     // HDG penalty
    PetscReal theta = .5;  // one-step theta method
    PetscReal T = 1, dt = 0, rtol = 1e-10, e_abs = 0, e_rel = 0;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, rnorm;
    const char* creason = NULL;
    PetscBool plot = true;
    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "timowave";
    char output_h5[PATH_MAX] = {0};
    char plot_scale[PATH_MAX] = "1";
    char domain_path[PATH_MAX] = "domains/single1.geo";
    PetscInt timowave_test = 0;
    PetscViewer h5_viewer = NULL;
    const char *pc_type;
    char name[64];

    (void)e_rel;

    PetscLogStage s_as, s_ts, s_rf;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    std::vector<PetscInt> itemp;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, errors, times;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsInt("-deg", "polynomial degree", NULL, poly_deg, &poly_deg, &is_set));
    PetscCall(PetscOptionsReal("-theta", "time-step averaging weight, 0 < theta <= 0.5, use theta=0.25 for CN", NULL, theta, &theta, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsInt("-nx", "number of refinements", NULL, nx, &nx, &is_set));
    PetscCall(PetscOptionsInt("-nt", "number of timesteps", NULL, nt, &nt, &is_set));
    PetscCall(PetscOptionsReal("-T", "end time", NULL, T, &T, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-plot", "plot solution", NULL, plot, &plot, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-sol_h5", "save solution as h5", NULL, output_h5, output_h5, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-domain", "domain path", NULL, domain_path, domain_path, PATH_MAX, &is_set));
    PetscCall(PetscOptionsInt("-test", "timowave test", NULL, timowave_test, &timowave_test, &is_set));
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

    dt = T / nt;

    HDGBase *hdg = NULL;
    PetscCall(PetscHDGCreate(poly_deg, timowave_test, domain_path, tau, theta, dt, &hdg));
    hdg->set_refinement(nx);

    PRIN2IY(timowave_test);
    PRIN2IY(poly_deg);
    PRIN2FY(tau);
    PRIN2FY(theta);
    PRIN2IY(nt);
    PRIN2IY(nx);
    PRIN2FY(dt);
    PRIN2FY(T);
    hdg->plot_option("fileName", output_filename);
    hdg->plot_option("outputDir", output_directory);
    hdg->plot_option("printFileNumber", "true");
    hdg->plot_option("scale", plot_scale);

    zero_v = hdg->zero_vector();
    N = zero_v.size();

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));
    PetscCall(VecSetBlockSize(rhs, hdg->n_dofs_per_node()));
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, nt+1, nt+1, &errors));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, nt+1, &times));
    PetscCall(PetscObjectSetName((PetscObject)times, "times"));

    temp = hdg->make_initial(zero_v);

    if (output_h5[0]) {
      PetscCall(PetscViewerHDF5Open(PETSC_COMM_SELF, output_h5, FILE_MODE_WRITE, &h5_viewer));
      PetscCall(PetscViewerHDF5PushGroup(h5_viewer, "/trace"));

      for (PetscInt k = 0; k < N; k++)
        PetscCall(VecSetValue(rhs, k, temp[k], INSERT_VALUES));
      PetscCall(VecAssemblyBegin(rhs));
      PetscCall(VecAssemblyEnd(rhs));
      snprintf(name, sizeof name, "timestep_%05d", 0);
      PetscCall(PetscObjectSetName((PetscObject)rhs, name));
      PetscCall(VecView(rhs, h5_viewer));
    }

    if (plot)
      hdg->plot_solution(temp, 0.);

    temp2 = hdg->errors(temp, 0);
    temp3 = hdg->norms(temp, 0);
    e_abs = PetscMax(temp2[0], e_abs);
    // e_rel = PetscMax(temp2[0] / temp3[0], e_rel);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e_abs0: %.5e\n", e_abs));
    PetscCall(VecSetValue(errors, 0, e_abs, INSERT_VALUES));
    PetscCall(VecSetValue(times,  0, 0, INSERT_VALUES));

    PRIN2S(s_as);
    mat_coo = hdg->trace_to_flux_mat(0.);
    PetscCall(MatCreateFromOptions(PETSC_COMM_SELF, NULL, 1, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(MatSetPreallocationCOO(mat, mat_coo.value_vec.size(), (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, PETSC_TRUE));
    PRIN2SP();

    PetscCall(PCRegister("net2as", PCCreate_Net2AS));

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(PCGetType(pc, &pc_type));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "pc_type: %s\n", pc_type));

    PRIN2S(s_ts);
    for (PetscInt i = 1; i <= nt; i++) {
      // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "------------ TIMESTEP %d -------\n", i));
      PetscReal ti = i*dt, error = 0;
        PetscCall(VecSetValue(times, i, ti, INSERT_VALUES));

        std::span<PetscReal> span;
        PetscCall(VecGetSpan(rhs, span));
        PetscLogStagePush(s_rf);
        hdg->residual_flux2(std::span{zero_v}, span, ti);
        PetscLogStagePop();

        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, rhs));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        hdg->set_data(span, ti);
        if (plot) hdg->plot_solution(span, ti);
        error = hdg->errors(span, ti)[0];
        e_abs = PetscMax(error, e_abs);
        PetscCall(VecSetValue(errors, i, error, INSERT_VALUES));

        if (h5_viewer) {
          snprintf(name, sizeof name, "timestep_%05d", (int)i);
          PetscCall(PetscObjectSetName((PetscObject)rhs, name));
          PetscCall(VecView(rhs, h5_viewer));
        }

        PetscCall(VecRestoreSpan(rhs, span));
    }
    PRIN2SP();

    PetscCall(VecAssemblyBegin(times));
    PetscCall(VecAssemblyEnd(times));

    if (h5_viewer) {
      PetscCall(PetscObjectSetName((PetscObject)rhs, "times"));
      PetscCall(VecView(times, h5_viewer));
      PetscCall(PetscViewerHDF5PopGroup(h5_viewer));
    }

    PetscCall(KSPGetConvergedReasonString(ksp, &creason));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));

    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));

    avg_iterations = ((PetscReal)iterations) / nt;

    PRIN2FY(e_abs);
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);
    PRIN2FY(avg_iterations);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "output: %s/%s.*.vtu\n", output_directory, output_filename));

    delete hdg;
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&times));
    PetscCall(VecDestroy(&errors));
    PetscCall(VecDestroy(&rhs));
    PetscCall(PetscViewerDestroy(&h5_viewer));

    PetscCall(PetscFinalize());
    return 0;
}
