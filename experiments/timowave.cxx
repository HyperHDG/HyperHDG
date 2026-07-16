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

template<unsigned int poly_deg, unsigned int stage, template<unsigned int, typename param_float_t> typename Test>
using HDGTimoWave = GlobalLoop::Hyperbolic<
  Topology::File<1,3>,
  Geometry::File<1,3>,
  NodeDescriptor::File<1,3>,
  LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2*poly_deg, Test, PetscReal, stage>,
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

//  instantiate the wave solver for a fixed test problem at runtime (degree, stages) choices.
// stages = 0 selects the matched coupling (deg 1/3/5 -> s = 1/2/3). The compiled combinations
// are kept sparse for compile-time reasons: the matched couplings plus the fixed-degree stage
// ladder (5, s) for temporal-order studies with the spatial discretization held constant
// (ne9-11). hdg must be deallocated with `delete`.
template<template<unsigned int, typename> typename Test>
static PetscErrorCode CreateDeg(
    PetscInt poly_deg, PetscInt stages, const char* path, PetscReal tau, PetscReal dt,
    HDGBase** hdg
) {
  PetscBool loc_lu_full = PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-loc_lu_full", &loc_lu_full, NULL));
  const std::vector<double> vals = {tau, dt, (double)loc_lu_full};
  PetscCall(InitTest<Test>(path));
  if (stages == 0)
    stages = (poly_deg + 1) / 2;  // matched: temporal order covers spatial order p+1
  switch (10 * poly_deg + stages) {
  case 11: *hdg = new HDGWrapper(HDGTimoWave<1,1,Test>(path, vals)); break;
  case 32: *hdg = new HDGWrapper(HDGTimoWave<3,2,Test>(path, vals)); break;
  case 51: *hdg = new HDGWrapper(HDGTimoWave<5,1,Test>(path, vals)); break;
  case 52: *hdg = new HDGWrapper(HDGTimoWave<5,2,Test>(path, vals)); break;
  case 53: *hdg = new HDGWrapper(HDGTimoWave<5,3,Test>(path, vals)); break;
  default:
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "no compiled instantiation for poly_deg = %d, stages = %d "
               "(available: (1,1), (3,2), (5,1), (5,2), (5,3))",
               (int)poly_deg, (int)stages);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// select the test problem by name. hdg must be deallocated with `delete`.
PetscErrorCode PetscHDGCreate(
    PetscInt poly_deg, PetscInt stages, const char* test,
    const char* path, PetscReal tau, PetscReal dt,
    HDGBase** hdg
) {
  PetscFunctionBeginUser;
  //if      (0 == strcmp(test, "stiffness")) PetscCall(CreateDeg<TimoshenkoStiffness>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "sinclamp")) PetscCall(CreateDeg<TimoshenkoSinClamp>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "gaussian")) PetscCall(CreateDeg<TimoshenkoGaussian>(poly_deg, stages, path, tau, dt, hdg));
  // commented out to cut compile time while the (deg, stages) ladder is compiled (ne9-11):
  // if (0 == strcmp(test, "drumhead")) PetscCall(CreateDeg<TimoshenkoDrumhead>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "constant")) PetscCall(CreateDeg<TimoshenkoConstant>(poly_deg, stages, path, tau, dt, hdg));
  if (0 == strcmp(test, "wave4"))     PetscCall(CreateDeg<TestTimoWave4>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave1"))    PetscCall(CreateDeg<TestTimoWave1>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave3"))    PetscCall(CreateDeg<TestTimoWave3>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "wave9"))    PetscCall(CreateDeg<TestTimoWave9>(poly_deg, stages, path, tau, dt, hdg));
  // else if (0 == strcmp(test, "clamped"))  PetscCall(CreateDeg<TimoWaveClamped>(poly_deg, stages, path, tau, dt, hdg));
  else PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "unknown test = \"%s\"", test);
  PetscFunctionReturn(PETSC_SUCCESS);
}

// all CLI-configurable run parameters
struct Options {
  PetscInt poly_deg = 1, stages = 0, nx = 1, nt = 1, tau_s = 0;
  PetscReal tau = 1;         // HDG penalty
  PetscReal T = 1, dt = 0;
  char plot[PATH_MAX] = {0};
  char plot_scale[PATH_MAX] = "1";
  char plot_values[64] = "all";
  char plot_props[64] = "all";
  PetscInt plot_stride = 1;
  PetscBool plot_energy = PETSC_TRUE;
  char domain_path[PATH_MAX] = "domains/single1.geo";
  char test[256] = {0};
  PetscBool print_timestep = PETSC_FALSE, mat_only = PETSC_FALSE, mem_max = PETSC_FALSE;
};

// parse CLI options; done is set when -help was requested and the run should stop
static PetscErrorCode ParseOptions(Options& o, PetscBool* done)
{
  PetscBool is_set;

  PetscFunctionBeginUser;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
  PetscCall(PetscOptionsInt("-deg", "polynomial degree", NULL, o.poly_deg, &o.poly_deg, &is_set));
  PetscCall(PetscOptionsInt("-stages", "number of Gauss stages (0 = matched to degree: 1/3/5 -> 1/2/3)", NULL, o.stages, &o.stages, &is_set));
  PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, o.tau, &o.tau, &is_set));
  PetscCall(PetscOptionsInt("-nx", "number of refinements", NULL, o.nx, &o.nx, &is_set));
  PetscCall(PetscOptionsInt("-nt", "number of timesteps", NULL, o.nt, &o.nt, &is_set));
  PetscCall(PetscOptionsReal("-T", "end time", NULL, o.T, &o.T, &is_set));
  PetscCall(PetscOptionsString("-plot", "plot solution using HyperHGD", NULL, o.plot, o.plot, PATH_MAX, &is_set));
  PetscCall(PetscOptionsInt("-plot_stride", "write plot fields every k-th timestep only (t=0 and the final step are always written); Steps/Values keeps the true times", NULL, o.plot_stride, &o.plot_stride, &is_set));
  PetscCall(PetscOptionsString("-plot_values", "PointData/values selection: all|disp|z|mag|none, or raw 'i,j,k' / 'mag:i,j,k' (Timoshenko layout: displacement = components 6,7,8)", NULL, o.plot_values, o.plot_values, sizeof(o.plot_values), &is_set));
  PetscCall(PetscOptionsString("-plot_props", "CellData/properties columns: all|beams|none or raw 'i,j,k'; beams = 7..15 = normals+widths+fiber_id (render with netvis --beams-cols 0,1,2:3,4,5:6:7 --beams-skip 8=-1)", NULL, o.plot_props, o.plot_props, sizeof(o.plot_props), &is_set));
  PetscCall(PetscOptionsBool("-plot_energy", "write per-cell CellData/energies at each plotted step", NULL, o.plot_energy, &o.plot_energy, &is_set));
  PetscCall(PetscOptionsString("-domain", "domain path", NULL, o.domain_path, o.domain_path, PATH_MAX, &is_set));
  PetscCall(PetscOptionsString("-test", "timowave test problem: wave4 (others commented out); -wave4_px sets its spatial phase", NULL, o.test, o.test, sizeof(o.test), &is_set));
  PetscCall(PetscOptionsBool("-print_timestep", "print timestep progress", NULL, o.print_timestep, &o.print_timestep, &is_set));
  PetscCall(PetscOptionsBool("-mat_only", "only compute matrix", NULL, o.mat_only, &o.mat_only, &is_set));
  PetscCall(PetscOptionsBool("-mem_max", "print memory stats in yaml", NULL, o.mem_max, &o.mem_max, &is_set));
  // dedicated flag: is_set is clobbered by every option parsed after -tau_s
  PetscBool tau_s_set = PETSC_FALSE;
  PetscCall(PetscOptionsInt("-tau_s", "hdg penalty parameter exponent s with tau ~ h^s", NULL, o.tau_s, &o.tau_s, &tau_s_set));
  if (tau_s_set) {
    switch (o.tau_s) {
    case  1: o.tau = 1./o.nx; break; // tau ~ h
    case  0: o.tau = 1;       break;
    case -1: o.tau = o.nx;    break; // tau ~ 1/h
    default:
      PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "tau_s expected 1,0,-1 found '%d'", o.tau_s);
    }
  }
  PetscOptionsEnd();

  o.dt = o.T / o.nt;

  PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", done, &is_set));
  if (*done)
    PetscCall(PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// log stages/events; file-static so helpers can log without them in every signature
static struct {
  PetscLogStage t2f, pa, ts, rf, mk, ksp;
  PetscLogEvent set, plot, errors;
} logs;

static PetscErrorCode RegisterLogStages()
{
  PetscFunctionBeginUser;
  PetscCall(PetscLogStageRegister("t2f", &logs.t2f));
  PetscCall(PetscLogStageRegister("ksp", &logs.ksp));
  PetscCall(PetscLogStageRegister("preallocation", &logs.pa));
  PetscCall(PetscLogStageRegister("make_initial", &logs.mk));
  PetscCall(PetscLogStageRegister("Timestepping", &logs.ts));
  PetscCall(PetscLogStageRegister("residual_flux", &logs.rf));
  PetscCall(PetscLogEventRegister("hdg_set", 0, &logs.set));
  PetscCall(PetscLogEventRegister("hdg_plot", 0, &logs.plot));
  PetscCall(PetscLogEventRegister("hdg_errors", 0, &logs.errors));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// YAML header identifying the run
static PetscErrorCode PrintRunHeader(const Options& o, PetscInt gauss_stages)
{
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "timowave_test: %s\n", o.test));
  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, "poly_deg", o.poly_deg));
  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, "tau", o.tau));
  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, "tau_s", o.tau_s));
  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, "nt", o.nt));
  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, "nx", o.nx));
  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, "dt", o.dt));
  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, "T", o.T));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "gauss_stages: %d\n", (int)gauss_stages));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// map the plot CLI presets onto HyperHDG plot options
static PetscErrorCode SetupPlotOptions(HDGBase* hdg, const Options& o)
{
  PetscFunctionBeginUser;
  hdg->plot_option("fileName", o.plot);
  hdg->plot_option("scale", o.plot_scale);
  hdg->plot_option("fileEnding", "vtkhdf");
  hdg->plot_option("energy", o.plot_energy ? "true" : "false");
  // presets for the 18-component Timoshenko values layout (netvis: disp = 6,7,8)
  if      (0 == strcmp(o.plot_values, "disp")) hdg->plot_option("valuesSelect", "6,7,8");
  else if (0 == strcmp(o.plot_values, "z"))    hdg->plot_option("valuesSelect", "8");
  else if (0 == strcmp(o.plot_values, "mag"))  hdg->plot_option("valuesSelect", "mag:6,7,8");
  else                                         hdg->plot_option("valuesSelect", o.plot_values);
  if      (0 == strcmp(o.plot_props, "beams")) hdg->plot_option("propertiesSelect", "7,8,9,10,11,12,13,14,15");
  else                                         hdg->plot_option("propertiesSelect", o.plot_props);
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Trace-space work vectors: distributed rhs/solution in the global numbering from
// distribute_domain, a per-rank local vector (owned + ghost) for the per-edge operations,
// and the scatter between them:
//  - SCATTER_FORWARD (INSERT): pull global dof values into sol_local for the per-edge operations
//    (make_initial / set_data / plot / errors), which read ghost endpoints.
//  - SCATTER_REVERSE (ADD): push locally-computed residual contributions to their owning rank.
struct TraceWork {
  Vec rhs = NULL, sol_local = NULL;
  VecScatter scatter = NULL;
  std::vector<PetscScalar> zero;  // zero trace input of the residual evaluations

  PetscErrorCode Create(HDGBase* hdg) {
    PetscInt bs = hdg->n_dofs_per_node();
    PetscInt N = hdg->size_of_system();
    PetscInt n_owned = hdg->n_owned_dofs();
    PetscInt n_local = hdg->n_local_dofs();

    PetscFunctionBeginUser;
    zero.assign(n_local, 0.);
    PetscCall(VecCreateFromOptions(PETSC_COMM_WORLD, NULL, bs, n_owned, N, &rhs));
    PetscCall(PetscObjectSetName((PetscObject)rhs, "trace"));

    auto gidx = hdg->local_to_global_dofs();   // global index of each local dof (length n_local)
    IS is_global;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_local, &sol_local));
    PetscCall(VecSetBlockSize(sol_local, bs));
    PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, n_local, (PetscInt*)gidx.data(),
                              PETSC_COPY_VALUES, &is_global));
    PetscCall(VecScatterCreate(rhs, is_global, sol_local, NULL, &scatter));
    PetscCall(ISDestroy(&is_global));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscErrorCode Destroy() {
    PetscFunctionBeginUser;
    PetscCall(VecDestroy(&rhs));
    PetscCall(VecDestroy(&sol_local));
    PetscCall(VecScatterDestroy(&scatter));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
};

// Combine a per-rank error (each rank L2-postprocesses over its owned edges: sqrt of a sum of
// squares) into the global error sqrt(sum_r local_r^2) via an all-reduce of the squared parts.
static std::vector<PetscReal> GlobalErrors(std::vector<PetscReal> e)
{
  for (auto& v : e) v *= v;
  MPI_Allreduce(MPI_IN_PLACE, e.data(), (PetscMPIInt)e.size(), MPIU_REAL, MPI_SUM,
                PETSC_COMM_WORLD);
  for (auto& v : e) v = PetscSqrtReal(v);
  return e;
}

// per-timestep error/norm history plus running maxima for the final report
struct ErrorTracker {
  Vec errors = NULL, norms = NULL;
  PetscReal e_abs = 0, n_abs = 0, e_trace = 0, n_trace = 0;

  PetscErrorCode Create(PetscInt nt) {
    PetscFunctionBeginUser;
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, nt+1, nt+1, &errors));
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "norm_", 1, nt+1, nt+1, &norms));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  // record global errors/norms at time t into slot step and the running maxima
  PetscErrorCode Record(HDGBase* hdg, std::span<PetscScalar> span, PetscReal t, PetscInt step) {
    PetscFunctionBeginUser;
    PetscLogEventBegin(logs.errors, 0,0,0,0);
    auto e = GlobalErrors(hdg->errors(span, t));
    PetscLogEventEnd(logs.errors, 0,0,0,0);
    auto n = GlobalErrors(hdg->norms(span, t));
    e_abs = PetscMax(e[0], e_abs);
    n_abs = PetscMax(n[0], n_abs);
    e_trace = PetscMax(e[1], e_trace);
    n_trace = PetscMax(n[1], n_trace);
    PetscCall(VecSetValue(errors, step, e[0], INSERT_VALUES));
    PetscCall(VecSetValue(norms, step, n[0], INSERT_VALUES));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscErrorCode Report() {
    PetscFunctionBeginUser;
    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));
    PetscCall(VecAssemblyBegin(norms));
    PetscCall(VecAssemblyEnd(norms));

    PetscReal e_rel = e_abs / n_abs;
    // relative error in the length-weighted skeleton norm (analytic trace norm from norms()[1]);
    // tests without an analytic solution (n_trace == 0) keep the absolute value
    if (n_trace > 0) e_trace /= n_trace;

    PRIN2FY(e_abs);
    PRIN2FY(n_abs);
    PRIN2FY(n_trace);
    PRIN2FY(e_rel);
    PRIN2FY(e_trace);
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscErrorCode Destroy() {
    PetscFunctionBeginUser;
    PetscCall(VecDestroy(&errors));
    PetscCall(VecDestroy(&norms));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
};

// Condensed stage operators, assembled once (time-constant), one per Gauss representative.
// At s = 1 this is exactly the classic trace operator. For s >= 2 one extra operator sits at
// the sentinel index n_gauss_reps(): the real SPD endpoint-algebraic system that recomputes
// the trace and the dual (n,m) from the advanced state after finalize_step (the value-form
// recombination only reaches temporal order ~s+1 for the algebraic variables).
static PetscErrorCode AssembleStageMats(HDGBase* hdg, std::vector<Mat>& mats)
{
  PetscInt bs = hdg->n_dofs_per_node();
  PetscInt N = hdg->size_of_system();
  PetscInt n_owned = hdg->n_owned_dofs();

  PetscFunctionBeginUser;
  mats.assign(hdg->n_gauss_reps() + (hdg->n_gauss_stages() > 1 ? 1 : 0), NULL);
  for (PetscInt rep = 0; rep < (PetscInt)mats.size(); rep++) {
    auto mat_coo = hdg->trace_to_flux_mat_stage(rep);
    mat_coo.eliminate_zeros();
    PetscInt ncoo = mat_coo.value_vec.size();
    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, n_owned, n_owned, N, N,
                                   &mats[rep]));
    PetscCall(MatSetPreallocationCOO(mats[rep], ncoo, (PetscInt*)mat_coo.row_vec.data(),
                                     (PetscInt*)mat_coo.col_vec.data()));
    // PETSc copied the COO indices; free them before staging the values (68 GB at net3)
    { auto drop_i = std::move(mat_coo.row_vec); }
    { auto drop_j = std::move(mat_coo.col_vec); }
    PetscCall(MatSetValuesCOO(mats[rep], (PetscScalar*)mat_coo.value_vec.data(),
                              INSERT_VALUES));

    // 3.24 retains the COO mapping for further MatSetValuesCOO calls that never come
    // (150-270 GB at net3) and has no API to drop it; swap into a clean duplicate instead
    {
      Mat mat_clean;
      PetscCall(MatDuplicate(mats[rep], MAT_COPY_VALUES, &mat_clean));
      PetscCall(MatDestroy(&mats[rep]));
      mats[rep] = mat_clean;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// One KSP per stage representative. s = 1 (real build) defaults to the classic CG + net2as;
// complex stage operators (s >= 2) are complex-symmetric, not Hermitian: direct LU by default.
static PetscErrorCode CreateStageKSPs(
    HDGBase* hdg, const std::vector<Mat>& mats, PetscReal rtol,
    KSPMonitorYAML_Ctx* monitor, std::vector<KSP>& ksps
) {
  PC pc;
  const char* pc_type;

  PetscFunctionBeginUser;
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PCRegister("net2as", PCCreate_Net2AS));
#endif
  PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

  ksps.assign(mats.size(), NULL);
  for (PetscInt rep = 0; rep < (PetscInt)mats.size(); rep++) {
    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, mats[rep], mats[rep]));
    PetscCall(KSPGetPC(ksp, &pc));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(PCSetType(pc, "net2as"));
#else
    PetscCall(KSPSetType(ksp, KSPPREONLY));
    PetscCall(PCSetType(pc, PCLU));
#endif
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    if (rep == 0)
      PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", monitor));
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
    ksps[rep] = ksp;
  }
  PetscCall(PetscTime(&monitor->t0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// One condensed solve at stage/endpoint index rep. Residual built from local parts: each rank
// evaluates over its owned edges into the local (owned + ghost) vector, then scatter-reverse
// with ADD assembles the distributed rhs; the solution is scattered forward so set_data sees
// full data on ghost endpoints.
static PetscErrorCode SolveStageSystem(
    HDGBase* hdg, KSP ksp, TraceWork& w,
    KSPMonitorYAML_Ctx* monitor, PetscInt step, PetscReal t, PetscInt rep, PetscInt* iterations
) {
  std::span<PetscScalar> span;
  PetscInt its;

  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(w.sol_local));
  PetscCall(VecGetSpan(w.sol_local, span));
  PetscLogStagePush(logs.rf);
  hdg->residual_flux_stage(std::span{w.zero}, span, rep, t);
  PetscLogStagePop();
  PetscCall(VecRestoreSpan(w.sol_local, span));

  PetscCall(VecZeroEntries(w.rhs));
  PetscCall(VecScatterBegin(w.scatter, w.sol_local, w.rhs, ADD_VALUES, SCATTER_REVERSE));
  PetscCall(VecScatterEnd(w.scatter, w.sol_local, w.rhs, ADD_VALUES, SCATTER_REVERSE));
  PetscCall(VecScale(w.rhs, -1.));

  // enorm monitoring (-ksp_monitor_yaml_enorm): the reference solve needs the assembled
  // rhs; only meaningful with -nt 1 (single-step conditioning studies, ne18-24)
  if (step == 1 && rep == 0)
    PetscCall(KSPMonitorYAML_Setup(ksp, w.rhs, monitor));

  PetscCall(KSPSolve(ksp, w.rhs, w.rhs));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  *iterations += its;

  PetscCall(VecScatterBegin(w.scatter, w.rhs, w.sol_local, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(w.scatter, w.rhs, w.sol_local, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGetSpan(w.sol_local, span));
  PetscLogEventBegin(logs.set, 0,0,0,0);
  hdg->set_data_stage(span, rep, t);  // stage trace zeta -> stash stage locals + trace
  PetscLogEventEnd(logs.set, 0,0,0,0);
  PetscCall(VecRestoreSpan(w.sol_local, span));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Advance state and endpoint trace by one step: one solve per stage representative (s = 1: the
// single real-valued stage; s >= 2: complex representatives, conjugate partners analytic), then
// finalize (state recombined per edge from the stashed stage solutions). For s >= 2 the
// algebraic variables (trace lambda, dual (n,m)) are then recomputed from the advanced state
// via the endpoint-algebraic solve at the sentinel index n_gauss_reps() -- the recombination
// extrapolates them at temporal order ~s+1 only, the static recovery restores order 2s. s = 1
// keeps the exact extrapolation (no extra solve on the production path).
static PetscErrorCode DoTimestep(
    HDGBase* hdg, const std::vector<KSP>& ksps, TraceWork& w,
    KSPMonitorYAML_Ctx* monitor, PetscInt step, PetscReal t, PetscInt* iterations
) {
  const PetscInt n_reps = hdg->n_gauss_reps();

  PetscFunctionBeginUser;
  for (PetscInt rep = 0; rep < n_reps; rep++)
    PetscCall(SolveStageSystem(hdg, ksps[rep], w, monitor, step, t, rep, iterations));

  hdg->finalize_step();

  if ((PetscInt)ksps.size() > n_reps)
    PetscCall(SolveStageSystem(hdg, ksps[n_reps], w, monitor, step, t, n_reps, iterations));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    PetscBool help = PETSC_FALSE;
    Options opt;
    PetscReal rtol = 1e-10;
    PetscInt iterations = 0;
    PetscReal avg_iterations = 0, rnorm = 0;
    const char* creason = NULL;
    KSPMonitorYAML_Ctx ksp_monitor_yaml_ctx;
    HDGBase* hdg = NULL;
    TraceWork work;
    ErrorTracker track;
    std::vector<Mat> stage_mats;
    std::vector<KSP> stage_ksps;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscCall(ParseOptions(opt, &help));
    if (help) {
      PetscFinalize();
      return 0;
    }

    if (opt.mem_max) PetscCall(PetscMemorySetGetMaximumUsage());
    PetscCall(RegisterLogStages());

    PetscCall(PetscHDGCreate(opt.poly_deg, opt.stages, opt.test, opt.domain_path, opt.tau, opt.dt, &hdg));
    // set_refinement rebuilds the hypernode factory and drops the distributed (owned/global) dof
    // numbering, so only refine when actually requested; nx==1 is the construction default and a
    // no-op for hyEdge_dim==1. Refinement is not yet supported together with the distributed
    // assembly path, so reject that combination instead of silently producing a wrong system.
    {
      PetscMPIInt comm_size;
      PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
      PetscCheck(opt.nx == 1 || comm_size == 1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
                 "refinement (-nx %" PetscInt_FMT ") is not supported on %d ranks; run with one rank",
                 opt.nx, comm_size);
      if (opt.nx != 1)
        hdg->set_refinement(opt.nx);
    }

#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(hdg->n_gauss_stages() == 1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "multi-stage gauss (deg 3) needs the complex-PETSc build: cmake preset 'complex'");
#endif

    PetscCall(PrintRunHeader(opt, hdg->n_gauss_stages()));
    PetscCall(SetupPlotOptions(hdg, opt));
    PetscCall(work.Create(hdg));
    PetscCall(track.Create(opt.nt));

    PRIN2S(logs.mk);
    {
      std::span<PetscScalar> span;
      hdg->make_initial(std::vector<PetscReal>(hdg->n_local_dofs(), 0.));
      PetscCall(VecGetSpan(work.sol_local, span));
      if (*opt.plot)
        hdg->plot_solution(span, 0.);
      PetscCall(track.Record(hdg, span, 0., 0));
      PetscCall(VecRestoreSpan(work.sol_local, span));
    }
    PRIN2SP();
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e_abs0: %.5e\n", track.e_abs));

    PRIN2S(logs.t2f);
    PetscCall(AssembleStageMats(hdg, stage_mats));
    PRIN2SP();

    PetscCall(MatPrintSymmetry("t2f_symmetry", stage_mats[0]));

    if (opt.mat_only) goto end;

    PRIN2S(logs.ksp);
    PetscCall(CreateStageKSPs(hdg, stage_mats, rtol, &ksp_monitor_yaml_ctx, stage_ksps));
    PRIN2SP();

    PRIN2S(logs.ts);
    for (PetscInt i = 1; i <= opt.nt; i++) {
      if (opt.print_timestep)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "#------------ TIMESTEP %d -------\n", i));
      PetscReal ti = i*opt.dt;

      PetscCall(DoTimestep(hdg, stage_ksps, work, &ksp_monitor_yaml_ctx, i, ti, &iterations));

      std::span<PetscScalar> span;
      PetscCall(VecGetSpan(work.sol_local, span));
      if (*opt.plot && (i % opt.plot_stride == 0 || i == opt.nt)) {
        PetscLogEventBegin(logs.plot, 0,0,0,0);
        hdg->plot_solution(span, ti);
        PetscLogEventEnd(logs.plot, 0,0,0,0);
      }
      PetscCall(track.Record(hdg, span, ti, i));
      PetscCall(VecRestoreSpan(work.sol_local, span));
    }
    PRIN2SP();

    if (!stage_ksps.empty() && stage_ksps[0]) {
      PetscCall(KSPGetConvergedReasonString(stage_ksps[0], &creason));
      PetscCall(KSPGetResidualNorm(stage_ksps[0], &rnorm));
    }

    avg_iterations = ((PetscReal)iterations) / opt.nt;

    PetscCall(track.Report());
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);
    PRIN2FY(avg_iterations);

end:
    PetscCall(PetscOptionsLeftYAML(NULL));

    delete hdg;
    PetscCall(track.Destroy());
    PetscCall(work.Destroy());
    for (auto& k : stage_ksps)
      PetscCall(KSPDestroy(&k));
    for (auto& m : stage_mats)
      PetscCall(MatDestroy(&m));

    if (opt.mem_max) {
      PetscLogDouble mem_max;
      PetscCall(PetscMemoryGetMaximumUsage(&mem_max));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mem_max: %.5e\n", mem_max));
    }

    PetscCall(PetscFinalize());
    return 0;
}
