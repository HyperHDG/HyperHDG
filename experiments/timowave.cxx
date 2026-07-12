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

// A test problem opts into runtime parameter setup by defining a static Init(path).
template<typename T>
concept HasInit = requires(const char* path) { T::Init(path); };

// Run a test problem's static-parameter setup if it provides Init(path); no-op otherwise.
template<template<unsigned int, typename> typename Test>
static PetscErrorCode InitTest(const char* path)
{
  PetscFunctionBeginUser;
  if constexpr (HasInit<Test<3, PetscReal>>)
    PetscCall(Test<3, PetscReal>::Init(path));
  PetscFunctionReturn(PETSC_SUCCESS);
}

//  instantiate the wave solver for a fixed test problem at a runtime polynomial degree.
// hdg must be deallocated with `delete`.
template<template<unsigned int, typename> typename Test>
static PetscErrorCode CreateDeg(
    PetscInt poly_deg, const char* path, PetscReal tau, PetscReal theta, PetscReal dt,
    HDGBase** hdg
) {
  PetscFunctionBeginUser;
  PetscCall(InitTest<Test>(path));
  switch (poly_deg) {
  case 1: *hdg = new HDGWrapper(HDGTimoWave<1,Test>(path, {tau, theta, dt})); break;
  case 2: *hdg = new HDGWrapper(HDGTimoWave<2,Test>(path, {tau, theta, dt})); break;
  case 3: *hdg = new HDGWrapper(HDGTimoWave<3,Test>(path, {tau, theta, dt})); break;
  case 6: *hdg = new HDGWrapper(HDGTimoWave<6,Test>(path, {tau, theta, dt})); break;
  default:
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "unsupported poly_deg = %d", (int)poly_deg);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// select the test problem by name. hdg must be deallocated with `delete`.
PetscErrorCode PetscHDGCreate(
    PetscInt poly_deg, const char* test,
    const char* path, PetscReal tau, PetscReal theta, PetscReal dt,
    HDGBase** hdg
) {
  PetscFunctionBeginUser;
  if      (0 == strcmp(test, "stiffness")) PetscCall(CreateDeg<TimoshenkoStiffness>(poly_deg, path, tau, theta, dt, hdg));
  else if (0 == strcmp(test, "sinclamp")) PetscCall(CreateDeg<TimoshenkoSinClamp>(poly_deg, path, tau, theta, dt, hdg));
  else if (0 == strcmp(test, "gaussian")) PetscCall(CreateDeg<TimoshenkoGaussian>(poly_deg, path, tau, theta, dt, hdg));
  else if (0 == strcmp(test, "drumhead")) PetscCall(CreateDeg<TimoshenkoDrumhead>(poly_deg, path, tau, theta, dt, hdg));
  else if (0 == strcmp(test, "wave4"))     PetscCall(CreateDeg<TestTimoWave4>(poly_deg, path, tau, theta, dt, hdg));
  else if (0 == strcmp(test, "constant")) PetscCall(CreateDeg<TimoshenkoConstant>(poly_deg, path, tau, theta, dt, hdg));
  else PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "unknown test = \"%s\"", test);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPrintSymmetry(const char* msg, Mat mat) {
  Mat AT, D;
  PetscReal nrm, nrm_a;

  PetscFunctionBeginUser;
  MatTranspose(mat, MAT_INITIAL_MATRIX, &AT);
  MatDuplicate(mat, MAT_COPY_VALUES, &D);
  MatAXPY(D, -1.0, AT, DIFFERENT_NONZERO_PATTERN);
  MatNorm(D, NORM_FROBENIUS, &nrm);
  MatNorm(mat, NORM_FROBENIUS, &nrm_a);
  PetscPrintf(PETSC_COMM_WORLD, "%s: %g\n", msg, (double)(nrm/nrm_a));
  MatDestroy(&AT); MatDestroy(&D);
  PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsLeftYAML(PetscOptions options) {
    PetscInt unused;
    char **names;
    char **values;

    PetscCall(PetscOptionsLeftGet(NULL, &unused, &names, &values));
    if (unused == 0) goto end;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# WARNING! There are options you set that were not used!\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "options_left:\n"));
    for (PetscInt i = 0; i < unused; i++)
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  - name: \"%s\"\n    value: \"%s\"\n", names[i], values[i]));
end:
    PetscCall(PetscOptionsLeftRestore(NULL, &unused, &names, &values));
    PetscCall(PetscOptionsSetValue(NULL, "-options_left", "0"));
    return 0;
}

int main(int argc, char **argv) {
    PetscBool help = false, is_set, print_timestep = PETSC_FALSE, set_mem_max = PETSC_FALSE;
    PetscInt nt = 1, nx = 1, poly_deg = 1, tau_s = 0;
    PetscInt N;            // global system size
    PetscReal tau = 1;     // HDG penalty
    PetscReal theta = .5;  // one-step theta method
    PetscReal T = 1, dt = 0, rtol = 1e-10, e_abs = 0, e_rel = 0, e_trace = 0;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, rnorm;
    const char* creason = NULL;
    PetscBool have_cache = PETSC_FALSE;
    char plot[PATH_MAX] = {0};
    char plot_scale[PATH_MAX] = "1";
    char plot_values[64] = "all";
    char plot_props[64] = "all";
    PetscInt plot_stride = 1;
    PetscBool plot_energy = PETSC_TRUE;
    char domain_path[PATH_MAX] = "domains/single1.geo";
    char mat_cache[PATH_MAX] = {0};
    char static_init[PATH_MAX] = {0};
    char timowave_test[256] = {0};
    const char *pc_type;
    KSPMonitorYAML_Ctx ksp_monitor_yaml_ctx;
    PetscBool mat_only = PETSC_FALSE;

    (void)e_rel;

    PetscLogStage s_t2f, s_pa, s_ts, s_rf, s_mk, s_ksp;
    PetscLogEvent e_set, e_plot, e_errors;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    std::vector<PetscInt> itemp;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs = NULL, sol_local = NULL, errors = NULL;
    VecScatter scatter = NULL;
    Mat mat;
    KSP ksp = NULL;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsInt("-deg", "polynomial degree", NULL, poly_deg, &poly_deg, &is_set));
    PetscCall(PetscOptionsReal("-theta", "one-step theta scheme, weight on the new time level; stiff modes amplify by -(1-theta)/theta per step, so theta >= 0.5: 0.5 = Crank-Nicolson (energy-conserving), 1 = fully implicit (damped); verified ne18-22", NULL, theta, &theta, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsInt("-nx", "number of refinements", NULL, nx, &nx, &is_set));
    PetscCall(PetscOptionsInt("-nt", "number of timesteps", NULL, nt, &nt, &is_set));
    PetscCall(PetscOptionsReal("-T", "end time", NULL, T, &T, &is_set));
    PetscCall(PetscOptionsString("-plot", "plot solution using HyperHGD", NULL, plot, plot, PATH_MAX, &is_set));
    PetscCall(PetscOptionsInt("-plot_stride", "write plot fields every k-th timestep only (t=0 and the final step are always written); Steps/Values keeps the true times", NULL, plot_stride, &plot_stride, &is_set));
    PetscCall(PetscOptionsString("-plot_values", "PointData/values selection: all|disp|z|mag|none, or raw 'i,j,k' / 'mag:i,j,k' (Timoshenko layout: displacement = components 6,7,8)", NULL, plot_values, plot_values, sizeof(plot_values), &is_set));
    PetscCall(PetscOptionsString("-plot_props", "CellData/properties columns: all|beams|none or raw 'i,j,k'; beams = 7..15 = normals+widths+fiber_id (render with netvis --beams-cols 0,1,2:3,4,5:6:7 --beams-skip 8=-1)", NULL, plot_props, plot_props, sizeof(plot_props), &is_set));
    PetscCall(PetscOptionsBool("-plot_energy", "write per-cell CellData/energies at each plotted step", NULL, plot_energy, &plot_energy, &is_set));
    PetscCall(PetscOptionsString("-mat_cache", "path to matrix cache", NULL, mat_cache, mat_cache, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-static", "path static init trace variables", NULL, static_init, static_init, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-domain", "domain path", NULL, domain_path, domain_path, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-test", "timowave test problem: stiffness, sinclamp, gaussian, drumhead, wave4, constant", NULL, timowave_test, timowave_test, sizeof(timowave_test), &is_set));
    PetscCall(PetscOptionsBool("-print_timestep", "print timestep progress", NULL, print_timestep, &print_timestep, &is_set));
    PetscCall(PetscOptionsBool("-mat_only", "only compute matrix", NULL, mat_only, &mat_only, &is_set));
    PetscCall(PetscOptionsInt("-tau_s", "set tau~h^s", NULL, tau_s, &tau_s, &is_set));
    PetscCall(PetscOptionsBool("-mem_max", "print memory stats in yaml", NULL, set_mem_max, &set_mem_max, &is_set));
    if (is_set) {
      PetscReal h = 1./nx;
      switch (tau_s) {
      case  1: tau = h;   break;
      case  0: tau = 1;   break;
      case -1: tau = 1/h; break;
      default:
        PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                   "expected tau_s = -1,0,1 found %lld", (long long)tau_s);
      }
    }
    PetscOptionsEnd();

    // theta < 0.5 amplifies every under-resolved stiff mode by |-(1-theta)/theta| > 1 per
    // step -- on the fiber networks (omega_max*dt >> 1) that is a guaranteed blow-up
    // (ne18-22: theta = 0.25 grew by exactly x3 per step). Warn, don't fail: on a problem
    // with resolved dynamics a conditionally stable theta may still be intentional.
    if (theta < 0.5)
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "# WARNING: theta = %g < 0.5 amplifies stiff modes by %g per step\n",
        (double)theta, (double)((1 - theta) / theta)));

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    if (set_mem_max) PetscCall(PetscMemorySetGetMaximumUsage());

    PetscCall(PetscLogStageRegister("t2f", &s_t2f));
    PetscCall(PetscLogStageRegister("ksp", &s_ksp));
    PetscCall(PetscLogStageRegister("preallocation", &s_pa));
    PetscCall(PetscLogStageRegister("make_initial", &s_mk));
    PetscCall(PetscLogStageRegister("Timestepping", &s_ts));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));
    PetscCall(PetscLogEventRegister("hdg_set", 0, &e_set));
    PetscCall(PetscLogEventRegister("hdg_plot", 0, &e_plot));
    PetscCall(PetscLogEventRegister("hdg_errors", 0, &e_errors));

    dt = T / nt;

    HDGBase *hdg = NULL;
    PetscCall(PetscHDGCreate(poly_deg, timowave_test, domain_path, tau, theta, dt, &hdg));
    // set_refinement rebuilds the hypernode factory and drops the distributed (owned/global) dof
    // numbering, so only refine when actually requested; nx==1 is the construction default and a
    // no-op for hyEdge_dim==1. Refinement is not yet supported together with the distributed
    // assembly path, so reject that combination instead of silently producing a wrong system.
    {
      PetscMPIInt comm_size;
      PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
      PetscCheck(nx == 1 || comm_size == 1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
                 "refinement (-nx %" PetscInt_FMT ") is not supported on %d ranks; run with one rank",
                 nx, comm_size);
      if (nx != 1)
        hdg->set_refinement(nx);
    }

    PRIN2SY(timowave_test);
    PRIN2IY(poly_deg);
    PRIN2FY(tau);
    PRIN2IY(tau_s);
    PRIN2FY(theta);
    PRIN2IY(nt);
    PRIN2IY(nx);
    PRIN2FY(dt);
    PRIN2FY(T);
    hdg->plot_option("fileName", plot);
    hdg->plot_option("scale", plot_scale);
    hdg->plot_option("fileEnding", "vtkhdf");
    hdg->plot_option("energy", plot_energy ? "true" : "false");
    // presets for the 18-component Timoshenko values layout (netvis: disp = 6,7,8)
    if      (0 == strcmp(plot_values, "disp")) hdg->plot_option("valuesSelect", "6,7,8");
    else if (0 == strcmp(plot_values, "z"))    hdg->plot_option("valuesSelect", "8");
    else if (0 == strcmp(plot_values, "mag"))  hdg->plot_option("valuesSelect", "mag:6,7,8");
    else                                       hdg->plot_option("valuesSelect", plot_values);
    if      (0 == strcmp(plot_props, "beams")) hdg->plot_option("propertiesSelect", "7,8,9,10,11,12,13,14,15");
    else                                       hdg->plot_option("propertiesSelect", plot_props);

    zero_v = hdg->zero_vector();              // local (owned + ghost), all zeros
    PetscInt bs = hdg->n_dofs_per_node();
    N = hdg->size_of_system();                // global system size
    PetscInt n_owned = hdg->n_owned_dofs();   // dofs owned by this rank
    PetscInt n_local = zero_v.size();         // n_local_dofs(): owned + ghost

    // Distributed system matrix (owned rows per rank) and matching distributed rhs. The KSP solves in
    // the global numbering from distribute_domain; per-edge local work happens in sol_local.
    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, n_owned, n_owned, N, N, &mat));
    PetscCall(MatCreateVecs(mat, NULL, &rhs));
    PetscCall(PetscObjectSetName((PetscObject)rhs, "trace"));

    // Per-rank local vector (owned + ghost) plus a scatter between it and the global rhs:
    //  - SCATTER_FORWARD (INSERT): pull global dof values into sol_local for the per-edge operations
    //    (make_initial_from_static / set_data / plot / errors), which read ghost endpoints.
    //  - SCATTER_REVERSE (ADD): push locally-computed residual contributions to their owning rank.
    {
      auto gidx = hdg->local_to_global_dofs();   // global index of each local dof (length n_local)
      IS is_global;
      PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_local, &sol_local));
      PetscCall(VecSetBlockSize(sol_local, bs));
      PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, n_local, (PetscInt*)gidx.data(),
                                PETSC_COPY_VALUES, &is_global));
      PetscCall(VecScatterCreate(rhs, is_global, sol_local, NULL, &scatter));
      PetscCall(ISDestroy(&is_global));
    }

    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, nt+1, nt+1, &errors));

    // Combine a per-rank error (each rank L2-postprocesses over its owned edges: sqrt of a sum of
    // squares) into the global error sqrt(sum_r local_r^2) via an all-reduce of the squared parts.
    auto global_errors = [&](std::vector<PetscReal> e) {
      for (auto& v : e) v *= v;
      MPI_Allreduce(MPI_IN_PLACE, e.data(), (PetscMPIInt)e.size(), MPIU_REAL, MPI_SUM,
                    PETSC_COMM_WORLD);
      for (auto& v : e) v = PetscSqrtReal(v);
      return e;
    };

    // Initial trace at t=0, built into the local vector sol_local (owned + ghost). For -static, the
    // file holds the global static trace written by `network` in the same global numbering: load it
    // into the distributed rhs, then scatter-forward so each rank has its full local part *including
    // ghost dofs* before make_initial_from_static reads per-edge (ghost-referencing) values.
    PRIN2S(s_mk);
    {
      std::span<PetscReal> span;
      if (*static_init) {
        PetscViewer viewer;
        PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, static_init, FILE_MODE_READ, &viewer));
        PetscCall(VecLoad(rhs, viewer));
        PetscCall(PetscViewerDestroy(&viewer));
        PetscCall(VecScatterBegin(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecGetSpan(sol_local, span));
        hdg->make_initial_from_static(span);
      }
      else {
        temp = hdg->make_initial(zero_v);
        PetscCall(VecGetSpan(sol_local, span));
        for (size_t k = 0; k < span.size(); ++k) span[k] = temp[k];
      }
      if (*plot)
        hdg->plot_solution(span, 0.);
      temp2 = global_errors(hdg->errors(span, 0));
      // temp3 = hdg->norms(...);  // per-rank only; would need an all-reduce for e_rel
      PetscCall(VecRestoreSpan(sol_local, span));
    }
    PRIN2SP();

    e_abs = PetscMax(temp2[0], e_abs);
    e_trace = PetscMax(temp2[1], e_trace);
    // e_rel = PetscMax(temp2[0] / temp3[0], e_rel);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e_abs0: %.5e\n", e_abs));
    PetscCall(VecSetValue(errors, 0, e_abs, INSERT_VALUES));

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
      // PETSc copies the index arrays into its own COO mapping, and MatSetValuesCOO only
      // needs the values -- free the indices here (~2/3 of the COO staging, 68 GB at net3)
      { auto drop_i = std::move(mat_coo.row_vec); }
      { auto drop_j = std::move(mat_coo.col_vec); }
      PetscCall(MatSetValuesCOO(mat, (PetscReal*)mat_coo.value_vec.data(), INSERT_VALUES));

      // PETSc retains internal COO mapping arrays on the matrix for repeated
      // MatSetValuesCOO calls that never come (the operator is time-constant), and
      // 3.24 has no API to drop them (~16-32 B per staged entry, 150-270 GB at net3).
      // Swap into a clean duplicate instead; transient cost is one extra matrix.
      {
        Mat mat_clean;
        PetscCall(MatDuplicate(mat, MAT_COPY_VALUES, &mat_clean));
        PetscCall(MatDestroy(&mat));
        mat = mat_clean;
      }
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


    PetscCall(MatPrintSymmetry("t2f_symmetry", mat));

    if (mat_only) goto end;

    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
    PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", &ksp_monitor_yaml_ctx));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(PCGetType(pc, &pc_type));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "pc_type: %s\n", pc_type));

    // Hand net2as the redistributed domain (points + edges in the partition's global numbering) so
    // its adjacency/points conform to the distributed system matrix instead of an independent file
    // read (which would use a different numbering).
    if (0 == strcmp(pc_type, "net2as")) {
      PetscInt node_bs = hdg->n_dofs_per_node();
      PetscInt n_owned_nodes = hdg->n_owned_dofs() / node_bs;
      PetscInt n_global_nodes = hdg->size_of_system() / node_bs;
      auto coords = hdg->owned_point_coords();
      auto edges_g = hdg->owned_edges_global();
      PetscCall(PCNet2ASSetDomain(pc, n_owned_nodes, n_global_nodes, hdg->n_space_dim(),
                                  coords.data(), (PetscInt)(edges_g.size() / 2),
                                  (PetscInt*)edges_g.data()));
    }

    PetscTime(&ksp_monitor_yaml_ctx.t0);
    PRIN2S(s_ksp);
    PetscCall(KSPSetUp(ksp));
    PRIN2SP();

    PRIN2S(s_ts);
    for (PetscInt i = 1; i <= nt; i++) {
      if (print_timestep)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "#------------ TIMESTEP %d -------\n", i));
      PetscReal ti = i*dt, error = 0;

        std::span<PetscReal> span;

        // Residual built from local parts: each rank evaluates the residual over its owned edges into
        // the local (owned + ghost) vector, then scatter-reverse with ADD pushes the ghost-row
        // contributions to their owning rank, assembling the distributed rhs.
        PetscCall(VecGetSpan(sol_local, span));
        PetscLogStagePush(s_rf);
        hdg->residual_flux2(std::span{zero_v}, span, ti);
        PetscLogStagePop();
        PetscCall(VecRestoreSpan(sol_local, span));

        PetscCall(VecZeroEntries(rhs));
        PetscCall(VecScatterBegin(scatter, sol_local, rhs, ADD_VALUES, SCATTER_REVERSE));
        PetscCall(VecScatterEnd(scatter, sol_local, rhs, ADD_VALUES, SCATTER_REVERSE));

        PetscCall(VecScale(rhs, -1.));

        // enorm monitoring (-ksp_monitor_yaml_enorm): the reference solve needs the
        // assembled rhs, so it happens here rather than next to KSPSetUp. The reference
        // is this first step's solution -- enorm is only meaningful with -nt 1
        // (single-step conditioning studies, ne18-24); later steps reuse a stale u_ref.
        if (i == 1)
          PetscCall(KSPMonitorYAML_Setup(ksp, rhs, &ksp_monitor_yaml_ctx));

        PetscCall(KSPSolve(ksp, rhs, rhs));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        // Pull the distributed solution into the local vector (owned + ghost) so the per-edge
        // operations below see full data on each rank's ghost endpoints.
        PetscCall(VecScatterBegin(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));

        PetscCall(VecGetSpan(sol_local, span));
        PetscLogEventBegin(e_set, 0,0,0,0);
        hdg->set_data(span, ti);
        PetscLogEventEnd(e_set, 0,0,0,0);
        if (*plot && (i % plot_stride == 0 || i == nt)) {
          PetscLogEventBegin(e_plot, 0,0,0,0);
          hdg->plot_solution(span, ti);
          PetscLogEventEnd(e_plot, 0,0,0,0);
        }
        PetscLogEventBegin(e_errors, 0,0,0,0);
        temp2 = global_errors(hdg->errors(span, ti));
        PetscLogEventEnd(e_errors, 0,0,0,0);
        error = temp2[0];
        e_abs = PetscMax(error, e_abs);
        e_trace = PetscMax(temp2[1], e_trace);
        PetscCall(VecRestoreSpan(sol_local, span));

        PetscCall(VecSetValue(errors, i, error, INSERT_VALUES));
    }
    PRIN2SP();

    PetscCall(KSPGetConvergedReasonString(ksp, &creason));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));

    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));

    avg_iterations = ((PetscReal)iterations) / nt;

    PRIN2FY(e_abs);
    PRIN2FY(e_trace);
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);
    PRIN2FY(avg_iterations);

end:
    PetscCall(PetscOptionsLeftYAML(NULL));

    delete hdg;
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&errors));
    PetscCall(VecDestroy(&rhs));
    PetscCall(VecDestroy(&sol_local));
    PetscCall(VecScatterDestroy(&scatter));

    if (set_mem_max) {
      PetscLogDouble mem_max;
      PetscCall(PetscMemoryGetMaximumUsage(&mem_max));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mem_max: %.5e\n", mem_max));
    }

    PetscCall(PetscFinalize());
    return 0;
}
