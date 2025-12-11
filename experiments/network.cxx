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
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec vertices;
  // bounding box of vertices
  PetscReal min[3], max[3];
  // number of subdomains in xy, total number of systems
  PetscInt p[2], sz;
  // coarse basis representation of the overlapping subdomains
  // rows correspond to subdomains
  Mat  sub;
  // mat[0] coarse system,
  // mat[i] for i >= 1 local matrices corresponding to subdomains
  Mat* mat;
  // corresponding solvers
  KSP* ksp;

  Vec sub_left;
};

PetscErrorCode PCDestroy_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscFunctionBegin;
  for (PetscInt i = 0; data->ksp && i < data->sz; i++) {
    PetscCall(KSPDestroy(data->ksp+i));
    PetscCall(MatDestroy(data->mat+i));
  }
  PetscCall(MatDestroy(&data->sub));
  PetscCall(VecDestroy(&data->sub_left));
  PetscCall(PetscFree(data->ksp));
  PetscCall(PetscFree(data->mat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool set;
  PetscInt p = data->p[0], p_lb = data->p[0];

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject, "Net2AS options");

  PetscCall(PetscOptionsBoundedInt("-pc_net2as_p", "number of subdomains per axis", NULL, p, &p, &set, p_lb));
  if (set) data->p[0] = data->p[1] = p;
  PetscCall(PetscOptionsBoundedInt("-pc_net2as_px", "number of subdomains", NULL, data->p[0], &data->p[0], &set, p_lb));
  PetscCall(PetscOptionsBoundedInt("-pc_net2as_py", "number of subdomains", NULL, data->p[1], &data->p[1], &set, p_lb));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec v = data->vertices;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt vstart, vend, size, nnz = 0, n_cols = data->p[0]*data->p[1], n_rows;
  PetscInt *rows, *cols;
  PetscReal *vals, h[2], eps = 1e-14;
  PetscBool done;
  const PetscInt *ioff, *inds;
  std::span<PetscReal> vspan;
  const char* prefix;
  PC subpc;
  Mat coarse_basis, A;
  MatType type;

  PetscFunctionBegin;
  PetscCall(PCDestroy_Net2AS(pc));
  PetscCall(PCGetOptionsPrefix(pc, &prefix));

  data->sz = n_cols+1;
  for (PetscInt i = 0; i < 3; i++) {
    PetscCall(VecStrideMin(v, i, NULL, data->min+i));
    PetscCall(VecStrideMax(v, i, NULL, data->max+i));
    if (i < 2) h[i] = (data->max[i]-data->min[i])/(data->p[i]+1);
  }

  PetscCall(VecGetOwnershipRange(v, &vstart, &vend));
  size = (vend-vstart)/3;
  PetscCall(PetscMalloc3(4zu*size, &rows, 4zu*size, &cols, 4zu*size, &vals));
  PetscCall(VecGetSpan(v, vspan));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "n -> x,y -> xx,yy -> i,j -> (type,col)\n"));
  for (PetscInt n = 0; n < size; n++) {
    PetscReal x = vspan[3*n],            y = vspan[3*n+1];
    PetscInt  i = (x-data->min[0])/h[0], j = (y-data->min[1])/h[1];
    // map to reference element
    PetscReal xx = (x-(i*h[0]+data->min[0]))/h[0], yy = (y-(j*h[1]+data->min[1]))/h[1];

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%3d -> %.3e,%.3e -> % .3e,% .3e -> % 3d,% 3d", n, x, y, xx, yy, i, j));

    if (i>data->p[0] || j>data->p[1]) goto next;

    if (i>0 && j>0) {
      rows[nnz] = vstart+n;
      cols[nnz] = (j-1)*data->p[0]+(i-1);
      vals[nnz] = (1-xx)*(1-yy);
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, " -> 0,%3d", cols[nnz]));
      nnz++;
    }

    if (i<data->p[0] && j>0) {
      rows[nnz] = vstart+n;
      cols[nnz] = (j-1)*data->p[0]+i;
      vals[nnz] = xx*(1-yy);
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, " -> 1,%3d", cols[nnz]));
      nnz++;
    }

    if (i>0 && j+1<data->p[1]) {
      rows[nnz] = vstart+n;
      cols[nnz] = j*data->p[0]+i-1;
      vals[nnz] = (1-xx)*yy;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, " -> 2,%3d", cols[nnz]));
      nnz++;
    }

    if (i+1<data->p[0] && j+1<data->p[1]) {
      rows[nnz] = vstart+n;
      cols[nnz] = j*data->p[0]+i;
      vals[nnz] = xx*yy;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, " -> 3,%3d", cols[nnz]));
      nnz++;
    }
  next:
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }
  PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "rows", rows, nnz));
  PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "cols", cols, nnz));
  PetscCall(PetscPrin2f(PETSC_COMM_WORLD, "vals", vals, nnz));

  PetscCall(VecRestoreSpan(v, vspan));
  PetscCall(PCGetOperators(pc, &A, NULL));
  PetscCall(MatGetType(A, &type));
  PetscCall(MatCreate(comm, &coarse_basis));
  PetscCall(MatSetType(coarse_basis, type));
  PetscCall(MatSetSizes(coarse_basis, size, n_cols, PETSC_DETERMINE, PETSC_DETERMINE));
  PetscCall(MatSetPreallocationCOO(coarse_basis, nnz, rows, cols));
  PetscCall(MatSetValuesCOO(coarse_basis, vals, INSERT_VALUES));
  PetscCall(MatEliminateZeros(coarse_basis, PETSC_TRUE));

  PetscCall(PetscMalloc2(data->sz, &data->ksp, data->sz, &data->mat));
  PetscCall(MatPtAP(A, coarse_basis, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->mat[0]));
  PetscCall(MatTranspose(coarse_basis, MAT_INITIAL_MATRIX, &data->sub)); // INITIAL -> redistributes
  PetscCall(MatView(data->sub, PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(MatFilter(data->sub, eps, /* compress = */ PETSC_TRUE, /* keep = */ PETSC_FALSE));
  PetscCall(MatView(data->sub, PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(MatCreateVecs(data->sub, NULL, &data->sub_left));

  PetscCall(MatGetRowIJ(data->sub, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");

  PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "ioff", ioff, n_rows+1));
  PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "inds", inds, ioff[n_rows]));

  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(KSPCreate(comm, data->ksp+i));

    PetscCall(KSPSetOptionsPrefix(data->ksp[i], prefix));
    PetscCall(KSPAppendOptionsPrefix(data->ksp[i], "net2as_"));

    PetscCall(KSPSetType(data->ksp[i], KSPPREONLY));
    PetscCall(KSPGetPC(data->ksp[i], &subpc));

    PetscCall(PCSetType(subpc, PCCHOLESKY));
    PetscCall(PCSetOptionsPrefix(subpc, prefix));
    PetscCall(PCAppendOptionsPrefix(subpc, "net2as_"));
    PetscCall(PCSetFromOptions(subpc));

    PetscCall(KSPSetFromOptions(data->ksp[i]));

    if (i > 0) {
      PetscInt s = i-1; // subdomain
      IS is;
      Mat *mat;
      PetscCall(ISCreateGeneral(PETSC_COMM_SELF, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_USE_POINTER, &is));
      // NOTE: MatCreateSubMatrix creates a submatrix of same type as A, regardless of comm of is,
      //       while MatCreateSubmatrices always creates sequential matrices,
      //       tough it also allocates the output parameter
      PetscCall(MatCreateSubMatrices(A, 1, &is, &is, MAT_INITIAL_MATRIX, &mat));
      data->mat[i] = *mat;
      PetscCall(PetscFree(mat));
      PetscCall(ISDestroy(&is));
    }
    PetscCall(KSPSetOperators(data->ksp[i], data->mat[i], data->mat[i]));
    PetscCall(KSPSetUp(data->ksp[i]));
  }
  PetscCall(MatRestoreRowIJ(data->sub, 0, PETSC_FALSE, PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  PetscCall(PetscFree3(rows, cols, vals));
  PetscCall(MatDestroy(&coarse_basis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));

  // coarse
  PetscCall(MatMult(data->sub, x, data->sub_left)); // sub = coarse_basis^T
  PetscCall(KSPSolve(data->ksp[0], data->sub_left, data->sub_left)); // override left
  PetscCall(MatMultTranspose(data->sub, data->sub_left, y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
  if (!isascii) goto end;

  PetscCall(PetscViewerASCIIPrintf(viewer, "p=%d,%d\n", data->p[0], data->p[1]));
  PetscCall(PetscViewerASCIIPrintf(viewer, "min=(%.5e,%.5e,%.5e)\n", data->min[0], data->min[1], data->min[2]));
  PetscCall(PetscViewerASCIIPrintf(viewer, "max=(%.5e,%.5e,%.5e)\n", data->max[0], data->max[1], data->max[2]));

  for (PetscInt i = 0; i < data->sz; i++) {
   PetscCall(PetscViewerASCIIPrintf(viewer, "--- SUB %d ---\n", i));
   PetscCall(MatView(data->mat[i], viewer));
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

  // minimal for Q1
  data->p[0] = data->p[1] = 2;

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
