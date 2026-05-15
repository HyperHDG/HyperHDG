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
  if (test == 0) {
    PetscViewer viewer;
    PetscReal size[3];
    PetscReal strain = .15;
    PetscInt  comp = 2;
    PetscBool is_set;

    PetscCall(PetscOptionsGetReal(NULL, NULL, "-strain", &strain, &is_set));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-comp", &comp, &is_set));
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5ReadAttribute(viewer, "/domain", "size", PETSC_DOUBLE, NULL, size));
    PetscCall(PetscViewerDestroy(&viewer));
    TimoshenkoStiffness<3>::length = size[0];
    TimoshenkoStiffness<3>::strain = strain;
    TimoshenkoStiffness<3>::comp = comp;
  }

  int i = poly_deg*10 + test;
  switch(i) {
  case 10: *hdg = new HDGWrapper(HDGTimoWave<1,TimoshenkoStiffness>(path, {tau, theta, dt})); return 0;
  case 14: *hdg = new HDGWrapper(HDGTimoWave<1,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 24: *hdg = new HDGWrapper(HDGTimoWave<2,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 34: *hdg = new HDGWrapper(HDGTimoWave<3,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  case 30: *hdg = new HDGWrapper(HDGTimoWave<3,TimoshenkoStiffness>(path, {tau, theta, dt})); return 0;
  case 64: *hdg = new HDGWrapper(HDGTimoWave<6,TestTimoWave4>(path, {tau, theta, dt})); return 0;
  default:
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "unsupported: poly_deg = %d, test = %d", poly_deg, test);
  }

  return 0;
}

int main(int argc, char **argv) {
    PetscBool help = false, is_set, print_timestep = PETSC_FALSE;
    PetscInt nt = 1, nx = 1, poly_deg = 1;
    PetscInt N;            // global system size
    PetscReal tau = 1;     // HDG penalty
    PetscReal theta = .5;  // one-step theta method
    PetscReal T = 1, dt = 0, rtol = 1e-10, e_abs = 0, e_rel = 0;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, rnorm;
    const char* creason = NULL;
    PetscBool have_cache = PETSC_FALSE;
    char plot[PATH_MAX] = {0};
    char plot_scale[PATH_MAX] = "1";
    char domain_path[PATH_MAX] = "domains/single1.geo";
    char mat_cache[PATH_MAX] = {0};
    char static_init[PATH_MAX] = {0};
    PetscInt timowave_test = 0;
    const char *pc_type;
    PetscBool ksp_monitor_yaml = PETSC_FALSE;
    KSPMonitorYAML_Ctx ksp_monitor_yaml_ctx;

    (void)e_rel;

    PetscLogStage s_t2f, s_pa, s_ts, s_rf, s_mk, s_ksp;

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
    PetscCall(PetscOptionsString("-plot", "plot solution using HyperHGD", NULL, plot, plot, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-mat_cache", "path to matrix cache", NULL, mat_cache, mat_cache, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-static", "path static init trace variables", NULL, static_init, static_init, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-domain", "domain path", NULL, domain_path, domain_path, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-ksp_monitor_yaml", "set yaml ksp monitor", NULL, ksp_monitor_yaml, &ksp_monitor_yaml, &is_set));
    PetscCall(PetscOptionsInt("-test", "timowave test", NULL, timowave_test, &timowave_test, &is_set));
    PetscCall(PetscOptionsBool("-print_timestep", "print timestep progress", NULL, print_timestep, &print_timestep, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    PetscCall(PetscLogStageRegister("t2f", &s_t2f));
    PetscCall(PetscLogStageRegister("ksp", &s_ksp));
    PetscCall(PetscLogStageRegister("preallocation", &s_pa));
    PetscCall(PetscLogStageRegister("make_initial", &s_mk));
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
    hdg->plot_option("fileName", plot);
    hdg->plot_option("scale", plot_scale);
    hdg->plot_option("fileEnding", "vtkhdf");

    zero_v = hdg->zero_vector();
    N = zero_v.size();

    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", hdg->n_dofs_per_node(), PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));
    PetscCall(VecSetBlockSize(rhs, hdg->n_dofs_per_node()));
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, nt+1, nt+1, &errors));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, nt+1, &times));
    PetscCall(PetscObjectSetName((PetscObject)times, "times"));
    PetscCall(PetscObjectSetName((PetscObject)rhs, "trace"));

    PRIN2S(s_mk);
    if (*static_init) {
      PetscViewer viewer;
      std::span<PetscReal> span;
      PetscCall(PetscViewerHDF5Open(PETSC_COMM_SELF, static_init, FILE_MODE_READ, &viewer));
      PetscCall(VecLoad(rhs, viewer));
      PetscCall(VecGetSpan(rhs, span));
      hdg->make_initial_from_static(span);
      temp.resize(span.size());
      std::copy(span.begin(), span.end(), temp.begin());
      PetscCall(VecRestoreSpan(rhs, span));
      PetscCall(VecDestroy(&rhs));
      PetscCall(PetscViewerDestroy(&viewer));
    }
    else {
      temp = hdg->make_initial(zero_v);
    }
    PRIN2SP();

    if (*plot)
      hdg->plot_solution(temp, 0.);

    return 1;

    temp2 = hdg->errors(temp, 0);
    temp3 = hdg->norms(temp, 0);
    e_abs = PetscMax(temp2[0], e_abs);
    // e_rel = PetscMax(temp2[0] / temp3[0], e_rel);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e_abs0: %.5e\n", e_abs));
    PetscCall(VecSetValue(errors, 0, e_abs, INSERT_VALUES));
    PetscCall(VecSetValue(times,  0, 0, INSERT_VALUES));

    PetscCall(PetscTestFile(mat_cache, 'r', &have_cache));
    if (have_cache) {
      PetscViewer viewer;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# loading matrix\n"));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mat_cache: %s\n", mat_cache));
      PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, mat_cache, FILE_MODE_READ, &viewer));
      PetscCall(MatLoad(mat, viewer));
      PetscCall(PetscViewerDestroy(&viewer));
    } else {
      PRIN2S(s_t2f);
      auto mat_coo = hdg->trace_to_flux_mat();
      mat_coo.eliminate_zeros();
      PetscInt ncoo = mat_coo.value_vec.size();
      PRIN2SP();

      PRIN2S(s_pa);
      PetscCall(MatSetPreallocationCOO(mat, ncoo, (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
      PetscCall(MatSetValuesCOO(mat, (PetscReal*)mat_coo.value_vec.data(), INSERT_VALUES));
      PRIN2SP();
    }

    if (!have_cache && *mat_cache) {
      PetscViewer viewer;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# saving matrix\n"));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mat_cache: %s\n", mat_cache));
      PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, mat_cache, FILE_MODE_WRITE, &viewer));
      PetscCall(MatView(mat, viewer));
      PetscCall(PetscViewerDestroy(&viewer));
    }

    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
    PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", &ksp_monitor_yaml_ctx));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(PCGetType(pc, &pc_type));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "pc_type: %s\n", pc_type));

    PetscTime(&ksp_monitor_yaml_ctx.t0);
    PRIN2S(s_ksp);
    PetscCall(KSPSetUp(ksp));
    PRIN2SP();

    PRIN2S(s_ts);
    for (PetscInt i = 1; i <= nt; i++) {
      if (print_timestep)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "------------ TIMESTEP %d -------\n", i));
      PetscReal ti = i*dt, error = 0;
        PetscCall(VecSetValue(times, i, ti, INSERT_VALUES));

        std::span<PetscReal> span;
        PetscCall(VecGetSpan(rhs, span));
        PetscLogStagePush(s_rf);
        hdg->residual_flux2(std::span{zero_v}, span, ti);
        PetscLogStagePop();
        PetscCall(VecRestoreSpan(rhs, span));

        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, rhs));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        PetscCall(VecGetSpan(rhs, span));
        hdg->set_data(span, ti);
        if (*plot) hdg->plot_solution(span, ti);
        error = hdg->errors(span, ti)[0];
        e_abs = PetscMax(error, e_abs);
        PetscCall(VecRestoreSpan(rhs, span));

        PetscCall(VecSetValue(errors, i, error, INSERT_VALUES));
    }
    PRIN2SP();

    PetscCall(VecAssemblyBegin(times));
    PetscCall(VecAssemblyEnd(times));

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
    if (*plot)
      PetscCall(PetscPrintf(PETSC_COMM_SELF, "output: output/%s.vtkhdf\n", plot));

    delete hdg;
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&times));
    PetscCall(VecDestroy(&errors));
    PetscCall(VecDestroy(&rhs));

    PetscCall(PetscFinalize());
    return 0;
}
