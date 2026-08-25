#include <petsc.h>
#include <stdio.h>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>
#include "hdg_base.hxx"
#include "net2as.hxx"
#include "parameters.hxx"
#include "prin2.hxx"

static const char help_msg[] = "experiments regarding timoshenko networks\n";

template <unsigned int deg, template <unsigned int, typename> typename Params>
using HDGElliptic = GlobalLoop::Elliptic<Topology::File<1, 3>,
                                         Geometry::File<1, 3>,
                                         NodeDescriptor::File<1, 3>,
                                         LocalSolver::TimoshenkoBeam<1, 3, deg, 2 * deg, Params> >;

// hdg must be deallocated with `delete`
PetscErrorCode PetscHDGCreate(const char* test, const char* domain, PetscReal tau, HDGBase** hdg)
{
  PetscFunctionBeginUser;

  if (0 == strcmp(test, "stiffness"))
  {
    PetscCall(TimoshenkoStiffness<>::Init(domain));
    *hdg = new HDGWrapper(HDGElliptic<3, TimoshenkoStiffness>(domain, tau));
    return 0;
  }
  else if (0 == strcmp(test, "gaussian"))
  {
    PetscCall(TimoshenkoGaussian<>::Init(domain));
    *hdg = new HDGWrapper(HDGElliptic<3, TimoshenkoGaussian>(domain, tau));
    return 0;
  }
  else
  {
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "unsupported test = '%s'", test);
  }

  return 0;
}

PetscErrorCode PetscOptionsLeftYAML(PetscOptions options)
{
  PetscInt unused;
  char** names;
  char** values;

  PetscCall(PetscOptionsLeftGet(NULL, &unused, &names, &values));
  if (unused == 0)
    goto end;

  PetscCall(
    PetscPrintf(PETSC_COMM_WORLD, "# WARNING! There are options you set that were not used!\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "options_left:\n"));
  for (PetscInt i = 0; i < unused; i++)
    PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "  - name: \"%s\"\n    value: \"%s\"\n", names[i], values[i]));
end:
  PetscCall(PetscOptionsLeftRestore(NULL, &unused, &names, &values));
  PetscCall(PetscOptionsSetValue(NULL, "-options_left", "0"));
  return 0;
}

PetscErrorCode MatPrintSymmetry(const char* msg, Mat mat)
{
  Mat AT, D;
  PetscReal nrm, nrm_a;

  PetscFunctionBeginUser;
  MatTranspose(mat, MAT_INITIAL_MATRIX, &AT);
  MatDuplicate(mat, MAT_COPY_VALUES, &D);
  MatAXPY(D, -1.0, AT, DIFFERENT_NONZERO_PATTERN);
  MatNorm(D, NORM_FROBENIUS, &nrm);
  MatNorm(mat, NORM_FROBENIUS, &nrm_a);
  PetscPrintf(PETSC_COMM_WORLD, "%s: %g\n", msg, (double)(nrm / nrm_a));
  MatDestroy(&AT);
  MatDestroy(&D);
  PetscFunctionReturn(0);
}

int main(int argc, char** argv)
{
  int rank, comm_size, proc_name_len;
  PetscReal rtol = 1e-10;

  char proc_name[MPI_MAX_PROCESSOR_NAME];
  char plot_path[PATH_MAX] = {0};
  char domain_filepath[PATH_MAX] = "domains/grid3.geo.bin";
  char viscoarse[PATH_MAX] = {0};
  char mat_cache[PATH_MAX] = {0};

  PetscLogStage s_as, s_it, s_rf, s_ksp, s_t2f, s_pa;

  PetscBool is_set, help, set_mem_max = PETSC_FALSE, mat_coo_off_proc = PETSC_FALSE;
  PetscInt N;
  PetscReal tau = 1;
  PetscInt iterations;
  PetscInt bs = 1;
  PetscReal emin, emax, cond;

  VecScatter scatter = NULL;
  Vec rhs = NULL, sol_local = NULL;
  Mat mat;
  KSP ksp;
  PC pc;
  PetscViewer viewer;
  const char* creason;
  PetscReal rnorm;
  PetscBool have_cache;
  PetscBool mat_only = PETSC_FALSE;
  PetscBool ksp_monitor_yaml = PETSC_FALSE;
  std::span<PetscReal> span;
  char test[10] = "stiffness";
  HDGBase* hdg = NULL;
  KSPMonitorYAML_Ctx ksp_monitor_yaml_ctx;
  const char* pc_type;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# ------- " __FILE__ " -------\n"));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
  PetscCallMPI(MPI_Get_processor_name(proc_name, &proc_name_len));
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Network Options", NULL);
  PetscCall(PetscOptionsString("-test", "test: (stiffness|gaussian)", NULL, test, test,
                               sizeof(test), &is_set));
  PetscCall(PetscOptionsString("-domain", "input network domain", NULL, domain_filepath,
                               domain_filepath, PATH_MAX, &is_set));
  PetscCall(PetscOptionsReal("-tau",
                             "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}",
                             NULL, tau, &tau, &is_set));
  PetscCall(PetscOptionsString("-plot", "plot solution using HyperHGD", NULL, plot_path, plot_path,
                               PATH_MAX, &is_set));
  PetscCall(PetscOptionsString("-viscoarse", "output name for visualization of coarse system", NULL,
                               viscoarse, viscoarse, PATH_MAX, &is_set));
  PetscCall(PetscOptionsBool("-mat_only", "only assemble matrix, overwrite any previous caches",
                             NULL, mat_only, &mat_only, &is_set));
  PetscCall(PetscOptionsString("-mat_cache", "path to matrix cache", NULL, mat_cache, mat_cache,
                               PATH_MAX, &is_set));
  PetscCall(PetscOptionsBool("-mem_max", "print memory stats in yaml", NULL, set_mem_max,
                             &set_mem_max, &is_set));
  PetscCall(PetscOptionsBool("-mat_coo_off_proc", "print memory stats in yaml", NULL,
                             mat_coo_off_proc, &mat_coo_off_proc, &is_set));
  PetscCall(PetscOptionsBool("-ksp_monitor_yaml", "set yaml ksp monitor", NULL, ksp_monitor_yaml,
                             &ksp_monitor_yaml, &is_set));
  PetscOptionsEnd();

  PetscCall(PetscPrin2Options());

  PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
  if (help)
  {
    Vec v;
    Mat m;
    KSP k;
    VecCreateFromOptions(PETSC_COMM_WORLD, NULL, 1, 1, 1, &v);
    MatCreateFromOptions(PETSC_COMM_WORLD, NULL, 1, 1, 1, 1, 1, &m);
    KSPCreate(PETSC_COMM_WORLD, &k);
    KSPSetFromOptions(k);
    VecDestroy(&v);
    MatDestroy(&m);
    KSPDestroy(&k);
    PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
    // PetscFinalize();
    // return 0;
  }

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mpi:\n  sz: %" PetscInt_FMT "\n", comm_size));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  names: ["));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s, ", proc_name));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "]\n"));

  if (set_mem_max)
    PetscCall(PetscMemorySetGetMaximumUsage());

  PetscCall(PetscLogStageRegister("assembly", &s_as));
  PetscCall(PetscLogStageRegister("iteration", &s_it));
  PetscCall(PetscLogStageRegister("residual", &s_rf));
  PetscCall(PetscLogStageRegister("ksp", &s_ksp));
  PetscCall(PetscLogStageRegister("t2f", &s_t2f));
  PetscCall(PetscLogStageRegister("prealloc", &s_pa));

  PetscCall(PetscHDGCreate(test, domain_filepath, tau, &hdg));
  PetscCall(PCRegister("net2as", PCCreate_Net2AS));
  PetscCall(
    KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetType(ksp, KSPCG));
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, "net2as"));
  PetscCall(KSPSetComputeSingularValues(ksp, PETSC_TRUE));
  PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
  PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", &ksp_monitor_yaml_ctx));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(PCGetType(pc, &pc_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "pc_type: %s\n", pc_type));

  bs = hdg->n_dofs_per_node();
  N = hdg->size_of_system();
  // Local row/col count = dofs owned by this rank (matches the partition's global numbering, so
  // PETSc's contiguous ownership ranges coincide with the renumbered owned dof blocks). Equals
  // PETSC_DECIDE behaviour for a single rank.
  PetscInt n_owned = hdg->n_owned_dofs();
  PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, n_owned, n_owned, N, N, &mat));
  PetscCall(KSPSetOperators(ksp, mat, mat));

  if (rank == 0)
    PetscCall(PetscTestFile(mat_cache, 'r', &have_cache));
  PetscCallMPI(MPI_Bcast(&have_cache, 1, MPI_C_BOOL, 0, PETSC_COMM_WORLD));
  if (have_cache)
  {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# loading matrix\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mat_cache: %s\n", mat_cache));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, mat_cache, FILE_MODE_READ, &viewer));
    PetscCall(MatLoad(mat, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  else
  {
    PRIN2S(s_t2f);
    auto mat_coo = hdg->trace_to_flux_mat();
    mat_coo.eliminate_zeros();
    PetscInt ncoo = mat_coo.value_vec.size();
    PRIN2SP();

    if (mat_coo_off_proc)
    {
      PetscInt off_count = 0, vstart, vend;
      PetscCall(MatGetOwnershipRange(mat, &vstart, &vend));
      for (PetscInt i = 0; i < ncoo; i++)
        if (mat_coo.row_vec[i] < (unsigned)vstart || mat_coo.row_vec[i] >= (unsigned)vend)
          off_count++;
      double off_frac = (double)off_count / ncoo;
      double off_max = 0;
      MPI_Reduce(&off_frac, &off_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mat_coo_off_proc: %.5e\n", off_max));
    }

    PRIN2S(s_pa);
    PetscCall(MatSetPreallocationCOO(mat, ncoo, (PetscInt*)mat_coo.row_vec.data(),
                                     (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, (PetscReal*)mat_coo.value_vec.data(), INSERT_VALUES));
    PRIN2SP();
  }

  if (!have_cache && *mat_cache)
  {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# saving matrix\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mat_cache: %s\n", mat_cache));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, mat_cache, FILE_MODE_WRITE, &viewer));
    PetscCall(MatView(mat, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }

  PetscCall(MatPrintSymmetry("t2f_symmetry", mat));

  {
    // Stage-1 validation: 1^T A 1 (= sum of all entries) and ||A 1|| are invariant under the
    // symmetric renumbering, so these must match between serial and distributed assembly.
    Vec ones, Aones;
    PetscReal ones_sum, ones_nrm;
    PetscCall(MatCreateVecs(mat, &ones, &Aones));
    PetscCall(VecSet(ones, 1.0));
    PetscCall(MatMult(mat, ones, Aones));
    PetscCall(VecSum(Aones, &ones_sum));
    PetscCall(VecNorm(Aones, NORM_2, &ones_nrm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t2f_ones_sum: %.12e\n", (double)ones_sum));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t2f_ones_nrm: %.12e\n", (double)ones_nrm));
    PetscCall(VecDestroy(&ones));
    PetscCall(VecDestroy(&Aones));
  }

  if (mat_only)
    goto end;
  PetscCall(MatCreateVecs(mat, NULL, &rhs));

  {
    // Map each local dof (owned then ghost) to its global index, and build a scatter that pulls
    // global dof values into a per-rank local vector (owned + ghost). Used to plot each rank's
    // owned edges (which reference ghost endpoints) from the distributed solution.
    auto gidx = hdg->local_to_global_dofs();
    PetscInt n_local = (PetscInt)gidx.size();
    IS is_global;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_local, &sol_local));
    PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, n_local, (PetscInt*)gidx.data(), PETSC_COPY_VALUES,
                              &is_global));
    PetscCall(VecScatterCreate(rhs, is_global, sol_local, NULL, &scatter));
    PetscCall(ISDestroy(&is_global));

    // Distributed residual: each rank computes the residual over its owned edges into a local
    // (owned + ghost) vector, then additively assembles into the global rhs by global dof index.
    // Off-process (ghost) rows are routed to their owners and summed by VecAssembly.
    PRIN2S(s_rf);
    auto zero = hdg->zero_vector();
    auto res_local = hdg->zero_vector();
    hdg->residual_flux2(zero, res_local, 0.);
    PetscCall(VecSetValues(rhs, n_local, (PetscInt*)gidx.data(), res_local.data(), ADD_VALUES));
    PetscCall(VecAssemblyBegin(rhs));
    PetscCall(VecAssemblyEnd(rhs));
    PRIN2SP();
  }

  PetscCall(VecScale(rhs, -1.));

  // Hand net2as the redistributed domain (points + edges in the partition's global numbering) so
  // its adjacency/points conform to the system matrix layout instead of an independent file read.
  if (0 == strcmp(pc_type, "net2as"))
  {
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

  PRIN2S(s_it);
  PetscCall(KSPSolve(ksp, rhs, rhs));
  PRIN2SP();

  PetscCall(KSPComputeExtremeSingularValues(ksp, &emax, &emin));
  cond = emax / emin;
  PetscCall(KSPGetIterationNumber(ksp, &iterations));
  PetscCall(KSPGetConvergedReasonString(ksp, &creason));
  PetscCall(KSPGetResidualNorm(ksp, &rnorm));
  PRIN2IY(iterations);
  PRIN2FY(rnorm);
  PRIN2FY(cond);
  PRIN2SY(creason);

  {
    // ||x|| is invariant under the symmetric renumbering, so it must match between serial and
    // distributed solves of the same physical problem.
    PetscReal solnorm;
    PetscCall(VecNorm(rhs, NORM_2, &solnorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "solnorm: %.12e\n", (double)solnorm));
  }

  if (*plot_path)
  {
    hdg->plot_option("fileName", plot_path);
    hdg->plot_option("fileNumber", "0");
    hdg->plot_option("fileEnding", "vtkhdf");
    PetscCall(VecScatterBegin(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, rhs, sol_local, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecGetSpan(sol_local, span));
    hdg->plot_solution(span);
    PetscCall(VecRestoreSpan(sol_local, span));
  }

  PetscCall(PetscObjectSetName((PetscObject)rhs, "trace"));
  PetscCall(VecViewFromOptions(rhs, NULL, "-trace_view"));

end:
  if (set_mem_max)
  {
    PetscLogDouble mem_max;
    PetscCall(PetscMemoryGetMaximumUsage(&mem_max));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mem_max: %.5e\n", mem_max));
  }

  PetscCall(PetscOptionsLeftYAML(NULL));

  PetscCall(KSPDestroy(&ksp));
  PetscCall(MatDestroy(&mat));
  PetscCall(VecDestroy(&rhs));
  PetscCall(VecDestroy(&sol_local));
  PetscCall(VecScatterDestroy(&scatter));

  PetscCall(PetscFinalize());
  return 0;
}
