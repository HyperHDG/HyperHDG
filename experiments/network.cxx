#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include "parameters.hxx"
#include "prin2.hxx"
#include "hdg_base.hxx"
#include "net2as.hxx"

static const char help_msg[] = "experiments regarding timoshenko networks\n";

template<unsigned int poly_deg>
using TB_LSol = LocalSolver::TimoshenkoBeam<1,3,poly_deg,2*poly_deg,LocalSolver::TimoschenkoBeamParametersClamped>;
template<unsigned int poly_deg>
using DF_LSol = LocalSolver::Diffusion<1,poly_deg,2*poly_deg,ConstantDiffusionParameters>;

template<unsigned int poly_deg, template<unsigned int> typename LSol>
using HDGNetwork = GlobalLoop::Elliptic<
  Topology::File<1,3>,
  Geometry::File<1,3>,
  NodeDescriptor::File<1,3>,
  LSol<poly_deg>
>;

// hdg must be deallocated with `delete`
PetscErrorCode PetscHDGCreate(const char* lsol, const char* domain, PetscReal tau, HDGBase **hdg) {
  if (0 == strcmp(lsol, "timo")) {
    *hdg = new HDGWrapper(HDGNetwork<3,TB_LSol>(domain, tau)); return 0;
  } else if (0 == strcmp(lsol, "diff")) {
    *hdg = new HDGWrapper(HDGNetwork<3,DF_LSol>(domain, tau)); return 0;
  } else {
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
      "unsupported lsol = %s", lsol);
  }

  return 0;
}

// NOTE: can only be called sequentially
PetscErrorCode PCNet2ASVisCoarse(PC pc, HDGBase* hdg, const char* name) {
  Mat cb;
  Vec left, right;
  PetscInt start, end;
  std::span<PetscReal> span;
  PCType type;

  PetscFunctionBeginUser;
  PCGetType(pc, &type);
  if (strcmp(type, "net2as") != 0) PetscFunctionReturn(0);
  PetscCall(PCNet2ASGetCB(pc, &cb));

  hdg->plot_option("fileName", name);
  hdg->plot_option("printFileNumber", "true");
  PetscCall(MatCreateVecs(cb, &right, &left));
  PetscCall(VecGetOwnershipRange(left, &start, &end));
  for (PetscInt i = 0; i < end-start; i++) {
    PetscCall(VecSetValue(left, i+start, 1., INSERT_VALUES));
    PetscCall(MatMultTranspose(cb, right, left));
    PetscCall(VecGetSpan(left, span));
    hdg->plot_solution(span, i);
    PetscCall(VecRestoreSpan(left, span));
    PetscCall(VecZeroEntries(right));
  }
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
    int rank, comm_size, proc_name_len;
    PetscReal rtol = 1e-10;

    char proc_name[MPI_MAX_PROCESSOR_NAME];
    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "network";
    char domain_filepath[PATH_MAX] = "domains/grid3.geo.bin";
    char plot_scale[PATH_MAX] = "1";
    char viscoarse[PATH_MAX] = {0};

    PetscLogStage s_as, s_it, s_rf, s_ksp;

    PetscBool is_set, help, set_mem_max = PETSC_FALSE;
    PetscInt N, ncoo;
    PetscReal tau = 1;
    PetscInt iterations;
    PetscInt bs = 1;

    sparse_mat<std::vector<PetscReal>> mat_coo;
    VecScatter scatter;
    Vec rhs, rhs0;
    Mat mat;
    KSP ksp;
    PC pc;
    PetscViewer viewer;
    const char* creason;
    PetscReal rnorm;
    PetscBool mat_only = PETSC_FALSE, mat_cached = PETSC_FALSE;
    PetscInt mat_load = 0;
    std::span<PetscReal> span;
    char lsol[10] = "timo";
    HDGBase* hdg = NULL;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# ------- " __FILE__ " -------\n"));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
    PetscCallMPI(MPI_Get_processor_name(proc_name, &proc_name_len));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Network Options", NULL);
    PetscCall(PetscOptionsString("-lsol", "local solver type: (timo|diff)", NULL, lsol, lsol, sizeof(lsol), &is_set));
    PetscCall(PetscOptionsString("-domain", "input network domain", NULL, domain_filepath, domain_filepath, PATH_MAX, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-viscoarse", "output name for visualization of coarse system", NULL, viscoarse, viscoarse, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-mat_only", "only assemble matrix, overwrite any previous", NULL, mat_only, &mat_only, &is_set));
    PetscCall(PetscOptionsInt("-mat_load", "force mat load if > 0, assembly if < 0", NULL, mat_load, &mat_load, &is_set));
    PetscCall(PetscOptionsBool("-mem_max", "print memory stats in yaml", NULL, set_mem_max, &set_mem_max, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscPrin2Options());

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      Vec v; Mat m; KSP k;
      VecCreateFromOptions(PETSC_COMM_WORLD, NULL, 1, 1, 1, &v);
      MatCreateFromOptions(PETSC_COMM_WORLD, NULL, 1, 1, 1, 1, 1, &m);
      KSPCreate(PETSC_COMM_WORLD, &k);
      KSPSetFromOptions(k);
      VecDestroy(&v);
      MatDestroy(&m);
      KSPDestroy(&k);
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mpi:\n  sz: %" PetscInt_FMT "\n", comm_size));
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  names:\n"));
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "    - %s\n", proc_name));

    if (set_mem_max) PetscCall(PetscMemorySetGetMaximumUsage());

    PetscCall(PetscLogStageRegister("assembly", &s_as));
    PetscCall(PetscLogStageRegister("iteration", &s_it));
    PetscCall(PetscLogStageRegister("residual", &s_rf));
    PetscCall(PetscLogStageRegister("ksp", &s_ksp));

    PetscCall(PetscHDGCreate(lsol, domain_filepath, tau, &hdg));
    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
    PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", NULL));
    PetscCall(KSPSetFromOptions(ksp));

    bs = hdg->n_dofs_per_node();
    N = hdg->size_of_system();
    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(KSPSetOperators(ksp, mat, mat));

    PRIN2S(s_as);
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, domain_filepath, FILE_MODE_APPEND, &viewer));
    PetscCall(PetscViewerHDF5HasGroup(viewer, "/mat", &mat_cached));
    PetscCall(PetscViewerHDF5PushGroup(viewer, "/mat"));
    if (mat_load == 0) mat_load = mat_cached ? 1 : -1;

    IS rows, cols;
    Vec vals;

    if (mat_load > 0) {
      PetscCall(ISCreate(PETSC_COMM_WORLD, &rows));
      PetscCall(ISCreate(PETSC_COMM_WORLD, &cols));
      PetscCall(VecCreate(PETSC_COMM_WORLD, &vals));
    } else {
      mat_coo = hdg->trace_to_flux_mat();
      ncoo = mat_coo.value_vec.size();
      PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, ncoo, (PetscInt*)mat_coo.row_vec.data(), PETSC_USE_POINTER, &rows));
      PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, ncoo, (PetscInt*)mat_coo.col_vec.data(), PETSC_USE_POINTER, &cols));
      PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, ncoo, PETSC_DETERMINE, mat_coo.value_vec.data(), &vals));
    }

    PetscCall(PetscObjectSetName((PetscObject)rows, "rows"));
    PetscCall(PetscObjectSetName((PetscObject)cols, "cols"));
    PetscCall(PetscObjectSetName((PetscObject)vals, "vals"));

    if (mat_load > 0) {
      const PetscReal *v;
      const PetscInt *r, *c;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "#   loading matrix\n"));

      PetscCall(ISLoad(rows, viewer));
      PetscCall(ISLoad(cols, viewer));
      PetscCall(VecLoad(vals, viewer));
      PetscCall(VecGetLocalSize(vals, &ncoo));
      mat_coo.resize(ncoo);

      PetscCall(ISGetIndices(rows, &r));
      PetscCall(ISGetIndices(cols, &c));
      PetscCall(VecGetArrayRead(vals, &v));

      PetscCall(PetscArraycpy((PetscInt*)mat_coo.row_vec.data(), r, ncoo));
      PetscCall(PetscArraycpy((PetscInt*)mat_coo.col_vec.data(), c, ncoo));
      PetscCall(PetscArraycpy((PetscReal*)mat_coo.value_vec.data(), v, ncoo));

      PetscCall(ISRestoreIndices(rows, &r));
      PetscCall(ISRestoreIndices(cols, &c));
      PetscCall(VecRestoreArrayRead(vals, &v));
    } else {
      PetscCall(ISView(rows, viewer));
      PetscCall(ISView(cols, viewer));
      PetscCall(VecView(vals, viewer));
    }

    // NOTE: may modify mat_coo
    PetscCall(MatSetPreallocationCOO(mat, ncoo, (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, (PetscReal*)mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, /* keep = */ PETSC_FALSE));

    PetscCall(ISDestroy(&rows));
    PetscCall(ISDestroy(&cols));
    PetscCall(VecDestroy(&vals));
    PetscCall(PetscViewerDestroy(&viewer));
    PRIN2SP();

    PetscCall(MatCreateVecs(mat, NULL, &rhs));

    if (mat_only) goto end;

    PRIN2S(s_ksp);
    PetscCall(KSPSetUp(ksp));
    PRIN2SP();

    if (strlen(viscoarse) > 0) PetscCall(PCNet2ASVisCoarse(pc, hdg, viscoarse));

    PetscCall(VecScatterCreateToZero(rhs, &scatter, &rhs0));
    if (rank == 0) {
      PRIN2S(s_rf);
      PetscCall(VecGetSpan(rhs0, span));
      auto zero = hdg->zero_vector();
      hdg->residual_flux2(zero, span, 0.);
      PetscCall(VecRestoreSpan(rhs0, span));
      PRIN2SP();
    }
    PetscCall(VecScatterBegin(scatter, rhs0, rhs, INSERT_VALUES, SCATTER_REVERSE));
    PetscCall(VecScatterEnd(scatter, rhs0, rhs, INSERT_VALUES, SCATTER_REVERSE));

    PetscCall(VecScale(rhs, -1.));

    PRIN2S(s_it);
    PetscCall(KSPSolve(ksp, rhs, rhs));
    PRIN2SP();

    PetscCall(KSPGetIterationNumber(ksp, &iterations));
    PetscCall(KSPGetConvergedReasonString(ksp, &creason));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);

    PetscCall(VecScatterBegin(scatter, rhs, rhs0, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, rhs, rhs0, INSERT_VALUES, SCATTER_FORWARD));

    hdg->plot_option("fileName", output_filename);
    hdg->plot_option("outputDir", output_directory);
    hdg->plot_option("printFileNumber", "false");
    hdg->plot_option("scale", plot_scale);
    if (rank == 0) {
      PetscCall(VecGetSpan(rhs0, span));
      hdg->plot_solution(span);
      PetscCall(VecRestoreSpan(rhs0, span));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "output: %s/%s.vtu\n", output_directory, output_filename));

    if (set_mem_max) {
      PetscLogDouble mem_max;
      PetscCall(PetscMemoryGetMaximumUsage(&mem_max));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mem_max: %.5e\n", mem_max));
    }

end:
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&rhs));
    PetscCall(VecDestroy(&rhs0));
    PetscCall(VecScatterDestroy(&scatter));

    PetscCall(PetscFinalize());
    return 0;
}
