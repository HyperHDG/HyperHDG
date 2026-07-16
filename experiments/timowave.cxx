#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>

#include <HyperHDG/local_solver/timowave.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>

#include <complex>

#include "parameters.hxx"
#include "hdg_base.hxx"
#include "petsc_util.hxx"
#if !defined(PETSC_USE_COMPLEX)
#include "net2as.hxx"  // real-build production preconditioner; not built against complex PETSc
#endif

static const char help_msg[] = "experiments regarding the wave equation\n";

// Gauss stage count tied to the spatial degree: temporal order 2s covers the spatial order
// p+1 with s = 1 for deg <= 2 and s = 2 (order 4, complex stage solves) for deg 3.
template<unsigned int poly_deg, template<unsigned int, typename param_float_t> typename Test>
using HDGTimoWave = GlobalLoop::Hyperbolic<
  Topology::File<1,3>,
  Geometry::File<1,3>,
  NodeDescriptor::File<1,3>,
  LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2*poly_deg, Test, PetscReal,
                              (poly_deg <= 2 ? 1u : 2u)>,
  std::vector<PetscScalar>
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
    PetscInt poly_deg, const char* path, PetscReal tau, PetscReal dt,
    HDGBase** hdg
) {
  PetscBool loc_lu_full = PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-loc_lu_full", &loc_lu_full, NULL));
  const std::vector<double> vals = {tau, dt, (double)loc_lu_full};
  PetscCall(InitTest<Test>(path));
  switch (poly_deg) {
  case 1: *hdg = new HDGWrapper(HDGTimoWave<1,Test>(path, vals)); break;
  case 2: *hdg = new HDGWrapper(HDGTimoWave<2,Test>(path, vals)); break;
  case 3: *hdg = new HDGWrapper(HDGTimoWave<3,Test>(path, vals)); break;
    //case 6: *hdg = new HDGWrapper(HDGTimoWave<6,Test>(path, vals)); break;
  default:
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "unsupported poly_deg = %d", (int)poly_deg);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// select the test problem by name. hdg must be deallocated with `delete`.
PetscErrorCode PetscHDGCreate(
    PetscInt poly_deg, const char* test,
    const char* path, PetscReal tau, PetscReal dt,
    HDGBase** hdg
) {
  PetscFunctionBeginUser;
  if      (0 == strcmp(test, "stiffness")) PetscCall(CreateDeg<TimoshenkoStiffness>(poly_deg, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "sinclamp")) PetscCall(CreateDeg<TimoshenkoSinClamp>(poly_deg, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "gaussian")) PetscCall(CreateDeg<TimoshenkoGaussian>(poly_deg, path, tau, dt, hdg));
  else if (0 == strcmp(test, "drumhead")) PetscCall(CreateDeg<TimoshenkoDrumhead>(poly_deg, path, tau, dt, hdg));
  else if (0 == strcmp(test, "wave4"))     PetscCall(CreateDeg<TestTimoWave4>(poly_deg, path, tau, dt, hdg));
  else if (0 == strcmp(test, "constant")) PetscCall(CreateDeg<TimoshenkoConstant>(poly_deg, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave1"))    PetscCall(CreateDeg<TestTimoWave1>(poly_deg, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave3"))    PetscCall(CreateDeg<TestTimoWave3>(poly_deg, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave9"))    PetscCall(CreateDeg<TestTimoWave9>(poly_deg, path, tau, dt, hdg));
  else if (0 == strcmp(test, "clamped"))  PetscCall(CreateDeg<TimoWaveClamped>(poly_deg, path, tau, dt, hdg));
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
    PetscReal T = 1, dt = 0, rtol = 1e-10, e_abs = 0, e_rel = 0, n_abs = 0, e_trace = 0,
              n_trace = 0;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, rnorm;
    const char* creason = NULL;
    char plot[PATH_MAX] = {0};
    char plot_scale[PATH_MAX] = "1";
    char plot_values[64] = "all";
    char plot_props[64] = "all";
    PetscInt plot_stride = 1;
    PetscBool plot_energy = PETSC_TRUE;
    char domain_path[PATH_MAX] = "domains/single1.geo";
    char timowave_test[256] = {0};
    const char *pc_type;
    KSPMonitorYAML_Ctx ksp_monitor_yaml_ctx;
    PetscBool mat_only = PETSC_FALSE;

    PetscLogStage s_t2f, s_pa, s_ts, s_rf, s_mk, s_ksp;
    PetscLogEvent e_set, e_plot, e_errors;

    std::vector<PetscReal> temp2, temp3;
    std::vector<PetscScalar> zero_s;  // zero trace input of the residual evaluations
    Vec rhs = NULL, sol_local = NULL, errors = NULL, norms = NULL;
    VecScatter scatter = NULL;
    // one condensed stage operator + KSP per representative (s = 1: the single real stage)
    std::vector<Mat> stage_mats;
    std::vector<KSP> stage_ksps;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsInt("-deg", "polynomial degree", NULL, poly_deg, &poly_deg, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsInt("-nx", "number of refinements", NULL, nx, &nx, &is_set));
    PetscCall(PetscOptionsInt("-nt", "number of timesteps", NULL, nt, &nt, &is_set));
    PetscCall(PetscOptionsReal("-T", "end time", NULL, T, &T, &is_set));
    PetscCall(PetscOptionsString("-plot", "plot solution using HyperHGD", NULL, plot, plot, PATH_MAX, &is_set));
    PetscCall(PetscOptionsInt("-plot_stride", "write plot fields every k-th timestep only (t=0 and the final step are always written); Steps/Values keeps the true times", NULL, plot_stride, &plot_stride, &is_set));
    PetscCall(PetscOptionsString("-plot_values", "PointData/values selection: all|disp|z|mag|none, or raw 'i,j,k' / 'mag:i,j,k' (Timoshenko layout: displacement = components 6,7,8)", NULL, plot_values, plot_values, sizeof(plot_values), &is_set));
    PetscCall(PetscOptionsString("-plot_props", "CellData/properties columns: all|beams|none or raw 'i,j,k'; beams = 7..15 = normals+widths+fiber_id (render with netvis --beams-cols 0,1,2:3,4,5:6:7 --beams-skip 8=-1)", NULL, plot_props, plot_props, sizeof(plot_props), &is_set));
    PetscCall(PetscOptionsBool("-plot_energy", "write per-cell CellData/energies at each plotted step", NULL, plot_energy, &plot_energy, &is_set));
    PetscCall(PetscOptionsString("-domain", "domain path", NULL, domain_path, domain_path, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-test", "timowave test problem: stiffness, sinclamp, gaussian, drumhead, wave4, constant", NULL, timowave_test, timowave_test, sizeof(timowave_test), &is_set));
    PetscCall(PetscOptionsBool("-print_timestep", "print timestep progress", NULL, print_timestep, &print_timestep, &is_set));
    PetscCall(PetscOptionsBool("-mat_only", "only compute matrix", NULL, mat_only, &mat_only, &is_set));
    PetscCall(PetscOptionsBool("-mem_max", "print memory stats in yaml", NULL, set_mem_max, &set_mem_max, &is_set));
    // dedicated flag: is_set is clobbered by every option parsed after -tau_s
    PetscBool tau_s_set = PETSC_FALSE;
    PetscCall(PetscOptionsInt("-tau_s", "hdg penalty parameter exponent s with tau ~ h^s", NULL, tau_s, &tau_s, &tau_s_set));
    if (tau_s_set) {
      switch (tau_s) {
      case  1: tau = 1./nx; break; // tau ~ h
      case  0: tau = 1;    break;
      case -1: tau = nx;   break; // tau ~ 1/h
      default:
        PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "tau_s expected 1,0,-1 found '%d'", tau_s);
      }
    }
    PetscOptionsEnd();

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
    PetscCall(PetscHDGCreate(poly_deg, timowave_test, domain_path, tau, dt, &hdg));
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

    // One unified time-stepping path for any stage count: per-representative stage operators
    // and solves. Multi-stage representatives are complex, so s >= 2 needs the complex-PETSc
    // build (cmake preset "complex"); s = 1 is real-valued and runs in both builds.
    const PetscInt n_greps = hdg->n_gauss_reps();
#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(hdg->n_gauss_stages() == 1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "multi-stage gauss (deg 3) needs the complex-PETSc build: cmake preset 'complex'");
#endif

    PRIN2SY(timowave_test);
    PRIN2IY(poly_deg);
    PRIN2FY(tau);
    PRIN2IY(tau_s);
    PRIN2IY(nt);
    PRIN2IY(nx);
    PRIN2FY(dt);
    PRIN2FY(T);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "gauss_stages: %d\n", (int)hdg->n_gauss_stages()));
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

    PetscInt bs = hdg->n_dofs_per_node();
    N = hdg->size_of_system();                // global system size
    PetscInt n_owned = hdg->n_owned_dofs();   // dofs owned by this rank
    PetscInt n_local = hdg->n_local_dofs();   // owned + ghost
    zero_s.assign(n_local, 0.);

    // Distributed rhs/solution vector in the global numbering from distribute_domain; per-edge
    // local work happens in sol_local.
    PetscCall(VecCreateFromOptions(PETSC_COMM_WORLD, NULL, bs, n_owned, N, &rhs));
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
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "norm_", 1, nt+1, nt+1, &norms));

    // Combine a per-rank error (each rank L2-postprocesses over its owned edges: sqrt of a sum of
    // squares) into the global error sqrt(sum_r local_r^2) via an all-reduce of the squared parts.
    auto global_errors = [&](std::vector<PetscReal> e) {
      for (auto& v : e) v *= v;
      MPI_Allreduce(MPI_IN_PLACE, e.data(), (PetscMPIInt)e.size(), MPIU_REAL, MPI_SUM,
                    PETSC_COMM_WORLD);
      for (auto& v : e) v = PetscSqrtReal(v);
      return e;
    };

    // Initial state at t=0: make_initial seeds the per-edge data (state coeffs_old AND the
    // endpoint trace lambda_old -- the trace is protocol-managed state, so the driver holds no
    // trace vector; diagnostics ignore their span argument). Static init (-static) was removed
    // pending its rework.
    PRIN2S(s_mk);
    {
      std::span<PetscScalar> span;
      hdg->make_initial(std::vector<PetscReal>(n_local, 0.));
      PetscCall(VecGetSpan(sol_local, span));
      if (*plot)
        hdg->plot_solution(span, 0.);
      temp2 = global_errors(hdg->errors(span, 0));
      temp3 = global_errors(hdg->norms(span, 0));
      PetscCall(VecRestoreSpan(sol_local, span));
    }
    PRIN2SP();

    e_abs = PetscMax(temp2[0], e_abs);
    n_abs = PetscMax(temp3[0], n_abs);
    e_trace = PetscMax(temp2[1], e_trace);
    n_trace = PetscMax(temp3[1], n_trace);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e_abs0: %.5e\n", e_abs));
    PetscCall(VecSetValue(errors, 0, e_abs, INSERT_VALUES));
    PetscCall(VecSetValue(norms,  0, temp3[0], INSERT_VALUES));

    // Condensed stage operators, assembled once (time-constant) and handed to one KSP per
    // representative. At s = 1 this is exactly the classic trace operator and the classic
    // KSP defaults (CG + net2as); complex stage systems (s >= 2) default to direct LU.
    PRIN2S(s_t2f);
    stage_mats.assign(n_greps, NULL);
    for (PetscInt rep = 0; rep < n_greps; rep++) {
      auto mat_coo = hdg->trace_to_flux_mat_stage(rep);
      mat_coo.eliminate_zeros();
      PetscInt ncoo = mat_coo.value_vec.size();
      PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, n_owned, n_owned, N, N,
                                     &stage_mats[rep]));
      PetscCall(MatSetPreallocationCOO(stage_mats[rep], ncoo, (PetscInt*)mat_coo.row_vec.data(),
                                       (PetscInt*)mat_coo.col_vec.data()));
      // PETSc copies the index arrays into its own COO mapping, and MatSetValuesCOO only
      // needs the values -- free the indices here (~2/3 of the COO staging, 68 GB at net3)
      { auto drop_i = std::move(mat_coo.row_vec); }
      { auto drop_j = std::move(mat_coo.col_vec); }
      PetscCall(MatSetValuesCOO(stage_mats[rep], (PetscScalar*)mat_coo.value_vec.data(),
                                INSERT_VALUES));

      // PETSc retains internal COO mapping arrays on the matrix for repeated
      // MatSetValuesCOO calls that never come (the operator is time-constant), and
      // 3.24 has no API to drop them (~16-32 B per staged entry, 150-270 GB at net3).
      // Swap into a clean duplicate instead; transient cost is one extra matrix.
      {
        Mat mat_clean;
        PetscCall(MatDuplicate(stage_mats[rep], MAT_COPY_VALUES, &mat_clean));
        PetscCall(MatDestroy(&stage_mats[rep]));
        stage_mats[rep] = mat_clean;
      }
    }
    PRIN2SP();

    PetscCall(MatPrintSymmetry("t2f_symmetry", stage_mats[0]));

    if (mat_only) goto end;

#if !defined(PETSC_USE_COMPLEX)
    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
#endif
    PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

    stage_ksps.assign(n_greps, NULL);
    PRIN2S(s_ksp);
    for (PetscInt rep = 0; rep < n_greps; rep++) {
      KSP ksp;
      PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
      PetscCall(KSPSetOperators(ksp, stage_mats[rep], stage_mats[rep]));
      PetscCall(KSPGetPC(ksp, &pc));
#if !defined(PETSC_USE_COMPLEX)
      PetscCall(KSPSetType(ksp, KSPCG));
      PetscCall(PCSetType(pc, "net2as"));
#else
      // the complex stage operators are complex-symmetric, not Hermitian: direct LU by default
      PetscCall(KSPSetType(ksp, KSPPREONLY));
      PetscCall(PCSetType(pc, PCLU));
#endif
      PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
      if (rep == 0)
        PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", &ksp_monitor_yaml_ctx));
      PetscCall(KSPSetFromOptions(ksp));
      if (rep == 0) {
        PetscCall(PCGetType(pc, &pc_type));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "pc_type: %s\n", pc_type));
      }

#if !defined(PETSC_USE_COMPLEX)
      // Hand net2as the redistributed domain (points + edges in the partition's global numbering)
      // so its adjacency/points conform to the distributed system matrix instead of an independent
      // file read (which would use a different numbering).
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
#endif
      PetscCall(KSPSetUp(ksp));
      stage_ksps[rep] = ksp;
    }
    PetscTime(&ksp_monitor_yaml_ctx.t0);
    PRIN2SP();

    PRIN2S(s_ts);
    for (PetscInt i = 1; i <= nt; i++) {
      if (print_timestep)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "#------------ TIMESTEP %d -------\n", i));
      PetscReal ti = i*dt, error = 0, norm = 0;

        std::span<PetscScalar> span;

        // One solve per stage representative (s = 1: the single real-valued stage; s >= 2:
        // complex representatives, conjugate partners analytic). Residual built from local
        // parts: each rank evaluates over its owned edges into the local (owned + ghost)
        // vector, then scatter-reverse with ADD assembles the distributed rhs; the solution is
        // scattered forward so set_data sees full data on ghost endpoints.
        for (PetscInt rep = 0; rep < n_greps; rep++) {
          PetscCall(VecZeroEntries(sol_local));
          PetscCall(VecGetSpan(sol_local, span));
          PetscLogStagePush(s_rf);
          hdg->residual_flux_stage(std::span{zero_s}, span, rep, ti);
          PetscLogStagePop();
          PetscCall(VecRestoreSpan(sol_local, span));

          PetscCall(VecZeroEntries(rhs));
          PetscCall(VecScatterBegin(scatter, sol_local, rhs, ADD_VALUES, SCATTER_REVERSE));
          PetscCall(VecScatterEnd(scatter, sol_local, rhs, ADD_VALUES, SCATTER_REVERSE));
          PetscCall(VecScale(rhs, -1.));

          // enorm monitoring (-ksp_monitor_yaml_enorm): the reference solve needs the assembled
          // rhs; only meaningful with -nt 1 (single-step conditioning studies, ne18-24)
          if (i == 1 && rep == 0)
            PetscCall(KSPMonitorYAML_Setup(stage_ksps[0], rhs, &ksp_monitor_yaml_ctx));

          PetscCall(KSPSolve(stage_ksps[rep], rhs, rhs));
          PetscCall(KSPGetIterationNumber(stage_ksps[rep], &its));
          iterations += its;

          PetscCall(VecScatterBegin(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
          PetscCall(VecScatterEnd(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
          PetscCall(VecGetSpan(sol_local, span));
          PetscLogEventBegin(e_set, 0,0,0,0);
          hdg->set_data_stage(span, rep, ti);  // stage trace zeta -> stash stage locals + trace
          PetscLogEventEnd(e_set, 0,0,0,0);
          PetscCall(VecRestoreSpan(sol_local, span));
        }

        // advances state AND endpoint trace (lambda+ = affine*lambda^n + sum mult*Re(w*zeta),
        // recombined per edge from the stashed stage traces)
        hdg->finalize_step();

        PetscCall(VecGetSpan(sol_local, span));
        if (*plot && (i % plot_stride == 0 || i == nt)) {
          PetscLogEventBegin(e_plot, 0,0,0,0);
          hdg->plot_solution(span, ti);
          PetscLogEventEnd(e_plot, 0,0,0,0);
        }
        PetscLogEventBegin(e_errors, 0,0,0,0);
        temp2 = global_errors(hdg->errors(span, ti));
        PetscLogEventEnd(e_errors, 0,0,0,0);
        error = temp2[0];
        temp3 = global_errors(hdg->norms(span, ti));
        norm = temp3[0];
        e_abs = PetscMax(error, e_abs);
        n_abs = PetscMax(norm, n_abs);
        e_trace = PetscMax(temp2[1], e_trace);
        n_trace = PetscMax(temp3[1], n_trace);
        PetscCall(VecRestoreSpan(sol_local, span));

        PetscCall(VecSetValue(errors, i, error, INSERT_VALUES));
        PetscCall(VecSetValue(norms, i, norm, INSERT_VALUES));
    }
    PRIN2SP();

    if (!stage_ksps.empty() && stage_ksps[0]) {
      PetscCall(KSPGetConvergedReasonString(stage_ksps[0], &creason));
      PetscCall(KSPGetResidualNorm(stage_ksps[0], &rnorm));
    }

    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));

    PetscCall(VecAssemblyBegin(norms));
    PetscCall(VecAssemblyEnd(norms));

    avg_iterations = ((PetscReal)iterations) / nt;
    e_rel = e_abs / n_abs;
    // relative error in the length-weighted skeleton norm (analytic trace norm from norms()[1]);
    // tests without an analytic solution (n_trace == 0) keep the absolute value
    e_trace = n_trace > 0 ? e_trace / n_trace : e_trace;

    PRIN2FY(e_abs);
    PRIN2FY(n_abs);
    PRIN2FY(n_trace);
    PRIN2FY(e_rel);
    PRIN2FY(e_trace);
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);
    PRIN2FY(avg_iterations);

end:
    PetscCall(PetscOptionsLeftYAML(NULL));

    delete hdg;
    PetscCall(VecDestroy(&errors));
    PetscCall(VecDestroy(&norms));
    PetscCall(VecDestroy(&rhs));
    PetscCall(VecDestroy(&sol_local));
    for (auto& k : stage_ksps)
      PetscCall(KSPDestroy(&k));
    for (auto& m : stage_mats)
      PetscCall(MatDestroy(&m));
    PetscCall(VecScatterDestroy(&scatter));

    if (set_mem_max) {
      PetscLogDouble mem_max;
      PetscCall(PetscMemoryGetMaximumUsage(&mem_max));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mem_max: %.5e\n", mem_max));
    }

    PetscCall(PetscFinalize());
    return 0;
}
