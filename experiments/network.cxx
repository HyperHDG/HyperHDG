#include <stdio.h>
#include <petsc.h>
#include <petsc/private/pcimpl.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include "parameters.hxx"
#include "geobin.hxx"

static const char help_msg[] = "experiments regarding timoshenko networks\n";
static PetscInt PETSC_PRIN2_ROW_LEN = 10;

PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, const PetscReal* dat, PetscInt len) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(com, msg));
  for (PetscInt i = 0; i < len; i++) {
    if (i % PETSC_PRIN2_ROW_LEN == 0)
       PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % .5e", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2i(MPI_Comm com, const char* msg, const PetscInt* dat, PetscInt len) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(com, msg));
  for (PetscInt i = 0; i < len; i++) {
    if (i % PETSC_PRIN2_ROW_LEN == 0)
      PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % 12d", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2Options() {
  PetscBool set;

  PetscFunctionBeginUser;
  PetscOptionsBegin(PETSC_COMM_WORLD, "prin2_", "Prin2", NULL);
  PetscCall(PetscOptionsInt("-row_len", "length of displayed rows", NULL, PETSC_PRIN2_ROW_LEN, &PETSC_PRIN2_ROW_LEN, &set));
  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

// must call VecRestoreSpan(x, span) after
PetscErrorCode VecGetSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p;
  PetscInt n;

  PetscFunctionBeginUser;
  PetscCall(VecGetArray(x, &p));
  PetscCall(VecGetLocalSize(x, &n));
  span = {p, (size_t)n};
  PetscFunctionReturn(0);
}

// must be called after each VecGetSpan(x, span)
PetscErrorCode VecRestoreSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p = span.data();

  PetscFunctionBeginUser;
  PetscCall(VecRestoreArray(x, &p));
  span = std::span<PetscScalar>();
  PetscFunctionReturn(0);
}

struct PC_Net2AS {
  Vec vertices;
  Mat coarse;
  PetscInt p;
  KSP* ksp;

  PetscReal min[3], max[3];
};

PetscErrorCode PCDestroy_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscFunctionBegin;
  PetscCall(MatDestroy(&data->coarse));
  for (PetscInt i = 0; data->ksp && i < data->p; i++)
    PetscCall(KSPDestroy(data->ksp+i));
  PetscFree(data->ksp);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool set;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject, "Net2AS options");
  PetscCall(PetscOptionsBoundedInt("-pc_net2as_p", "number of subdomains", NULL, data->p, &data->p, &set, data->p));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec v = data->vertices;
  PetscInt vstart, vend, size;
  PetscInt *cooi, *cooj;
  PetscReal *coov;
  std::span<PetscReal> vspan;
  const char* prefix;
  PC subpc;

  PetscFunctionBegin;
  PetscCall(PCDestroy_Net2AS(pc));
  PetscCall(PCGetOptionsPrefix(pc, &prefix));

  PetscCall(VecGetSpan(v, vspan));
  PetscCall(VecGetOwnershipRange(v, &vstart, &vend));
  size = vend-vstart;
  PetscCall(PetscMalloc3(4*size, &cooi, 4*size, &cooj, 4*size, &coov));
  for (PetscInt i = 0; i < 3; i++) {
    PetscCall(VecStrideMin(v, i, NULL, data->min+i));
    PetscCall(VecStrideMax(v, i, NULL, data->max+i));
  }
  PetscCall(VecRestoreSpan(data->vertices, vspan));

  PetscCall(PetscMalloc1(data->p+1, &data->ksp));
  for (PetscInt i = 0; i < data->p+1; i++) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)pc), data->ksp+i));

    PetscCall(KSPSetOptionsPrefix(data->ksp[i], prefix));
    PetscCall(KSPAppendOptionsPrefix(data->ksp[i], "net2as_"));

    PetscCall(KSPSetType(data->ksp[i], KSPPREONLY));
    PetscCall(KSPGetPC(data->ksp[i], &subpc));

    PetscCall(PCSetType(subpc, PCCHOLESKY));
    PetscCall(PCSetOptionsPrefix(subpc, prefix));
    PetscCall(PCAppendOptionsPrefix(subpc, "net2as_"));
    PetscCall(PCSetFromOptions(subpc));

    PetscCall(KSPSetFromOptions(data->ksp[i]));
    // PetscCall(KSPSetUp(data->ksp[i]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));
  PetscCall(VecCopy(x, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
  if (!isascii) goto end;

  PetscCall(PetscViewerASCIIPrintf(viewer, "p=%d\n", data->p));
  PetscCall(PetscViewerASCIIPrintf(viewer, "min=(%.5e,%.5e,%.5e)\n", data->min[0], data->min[1], data->min[2]));
  PetscCall(PetscViewerASCIIPrintf(viewer, "max=(%.5e,%.5e,%.5e)\n", data->max[0], data->max[1], data->max[2]));

  for (PetscInt i = 0; i < data->p+1; i++) {
   PetscCall(PetscViewerASCIIPrintf(viewer, "sub KSP %d\n", i));
   PetscCall(KSPView(data->ksp[i], viewer));
  }
end:
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCNet2ASReadGraph(PC pc, const char *path) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  std::span<PetscReal> vspan;
  PetscInt start, end;
  PetscReal *verts;

  PetscFunctionBeginUser;
  geobin::Graph graph = geobin::deserialize_bin(path);
  verts = (PetscReal*)graph.vertices.data();
  PetscCall(VecCreateFromOptions(PetscObjectComm((PetscObject)pc), NULL, 3, PETSC_DECIDE, graph.vertices.size()*3, &data->vertices));
  PetscCall(VecGetSpan(data->vertices, vspan));
  PetscCall(VecGetOwnershipRange(data->vertices, &start, &end));
  std::copy(verts+start, verts+end, vspan.data());
  PetscCall(VecRestoreSpan(data->vertices, vspan));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCreate_Net2AS(PC pc) {
  PC_Net2AS *data;

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&data));
  pc->data = (void*)data;

  data->p = 2; // minimal for Q1

  pc->ops->apply = PCApply_Net2AS;
  pc->ops->setup = PCSetup_Net2AS;
  pc->ops->destroy = PCDestroy_Net2AS;
  pc->ops->setfromoptions = PCSetFromOptions_Net2AS;
  pc->ops->view = PCView_Net2AS;

  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    // constexpr unsigned int poly_deg = 5;
    using Top = Topology::File<1,3>;
    using Geo = Geometry::File<1,3>;
    using NDes = NodeDescriptor::File<1,3>;
    // using LSol = LocalSolver::TimoshenkoBeam<1,3,poly_deg,2*poly_deg,LocalSolver::TimoschenkoBeamParametersClamped>;
    using LSol = LocalSolver::Diffusion<1,5,10,ConstantDiffusionParameters>;
    using HDG = GlobalLoop::Elliptic<Top,Geo,NDes,LSol>;
    constexpr PetscInt n_dofs_per_node = LSol::n_glob_dofs_per_node();

    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "network";
    char domain_filepath[PATH_MAX] = "domains/grid_8.geo.bin";
    char plot_scale[PATH_MAX] = "1";

    PetscLogStage s_as, s_it, s_rf;

    PetscBool is_set, help;
    PetscInt N, ncoo;
    PetscReal tau = 1;
    PetscInt iterations;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Network Options", NULL);
    PetscCall(PetscOptionsString("-domain", "input network domain", NULL, domain_filepath, domain_filepath, PATH_MAX, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscPrin2Options());

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "n_dofs_per_node", &n_dofs_per_node, 1));

    PetscCall(PetscLogStageRegister("Assembly", &s_as));
    PetscCall(PetscLogStageRegister("Iteration", &s_it));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));

    HDG hdg(domain_filepath);
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "false");
    hdg.plot_option("scale", plot_scale);

    zero_v = hdg.zero_vector();
    N = zero_v.size();
    PetscCall(VecCreate(PETSC_COMM_WORLD, &sol));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &rhs));
    PetscCall(VecSetType(sol, VECMPI));
    PetscCall(VecSetType(rhs, VECMPI));
    PetscCall(VecSetSizes(sol, PETSC_DECIDE, N));
    PetscCall(VecSetSizes(rhs, PETSC_DECIDE, N));

    PetscLogStagePush(s_as);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "assembly...\n"));
    mat_coo = hdg.trace_to_flux_mat();
    ncoo = mat_coo.value_vec.size();
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
    PetscCall(MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, N, N));
    PetscCall(MatSetType(mat, MATMPIAIJ));
    PetscCall(MatSetPreallocationCOO(mat, ncoo, (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
    PetscLogStagePop();

    PetscCall(PCRegister("net2as", PCCreate_Net2AS));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(PCNet2ASReadGraph(pc, domain_filepath));
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iteration...\n"));
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;
    PetscCall(VecGetSpan(rhs, rhs_span));
    PetscCall(VecGetSpan(sol, sol_span));

    PetscLogStagePush(s_rf);
    hdg.residual_flux2(zero_v, rhs_span, 0.);
    PetscLogStagePop();
    PetscCall(VecScale(rhs, -1.));

    PetscLogStagePush(s_it);
    PetscCall(KSPSolve(ksp, rhs, sol));
    PetscLogStagePop();

    PetscCall(KSPGetIterationNumber(ksp, &iterations));

    hdg.plot_solution(sol_span);

    PetscCall(VecRestoreSpan(rhs, rhs_span));
    PetscCall(VecRestoreSpan(sol, sol_span));

    PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "iterations", &iterations, 1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "wrote output to '%s/%s.*.vtu'\n", output_directory, output_filename));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}
