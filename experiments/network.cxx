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
static PetscLogDouble PETSC_PRIN2_TIMER = 0;
static PetscInt PETSC_PRIN2_STAGE = 0;
static const char* PETSC_PRIN2_STAGE_NAME = "";

#define PRIN2IY(VAR)  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2FY(VAR)  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2S(STAGE) do { PETSC_PRIN2_STAGE = STAGE; PetscCall(PetscLogStageGetName(STAGE, &PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# %s...\n", PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscTime(&PETSC_PRIN2_TIMER)); PetscCall(PetscLogStagePush(STAGE)); } while(0)
#define PRIN2SP()     do { PetscLogDouble time; PetscCall(PetscLogStagePop()); PetscCall(PetscTime(&time)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t_%s: %.5e\n", PETSC_PRIN2_STAGE_NAME, (time-PETSC_PRIN2_TIMER))); } while(0)

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

PetscErrorCode PetscPrin2iy(MPI_Comm comm, const char *name, PetscInt val) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: %d\n", name, val));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2fy(MPI_Comm comm, const char *name, PetscReal val) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: %.5e\n", name, val));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2iya(MPI_Comm comm, const char *name, PetscInt val, PetscInt num) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: [", name));
  for (PetscInt i = 0; i < num; i++)
    PetscCall(PetscPrintf(comm, "%d, ", val));
  PetscCall(PetscPrintf(comm, "]\n"));
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

PetscErrorCode KSPMonitorCSVCreate(PetscViewer viewer, PetscViewerFormat format, void *ctx, PetscViewerAndFormat **vf) {
    PetscFunctionBegin;
    PetscCall(PetscViewerAndFormatCreate(viewer, format, vf));
    PetscCall(PetscViewerASCIIPrintf(viewer, "iteration,residual_norm\n"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPMonitorCSVDestroy(PetscViewerAndFormat **vf) {
    PetscFunctionBegin;
    PetscCall(PetscViewerAndFormatDestroy(vf));
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode KSPMonitorCSV(KSP ksp, PetscInt it, PetscReal rnorm, PetscViewerAndFormat *vf) {
  PetscViewer viewer = vf->viewer;

  PetscFunctionBegin;
  PetscCall(PetscViewerASCIIPrintf(viewer, "%03" PetscInt_FMT ",%.16e\n", it, (double)rnorm));
  PetscFunctionReturn(0);
}

struct PC_Net2AS {
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec vertices;
  // bounding box of vertices
  PetscReal min[3], max[3];
  // number of subdomains in xy, total number of systems, block size (number of dofs per node)
  PetscInt p[2], sz, bs;
  // coarse basis representation of the overlapping subdomains
  // rows correspond to subdomains
  Mat  sub;
  // mat[0] coarse system,
  // mat[i] for i >= 1 local matrices corresponding to subdomains
  Mat* mat;
  // corresponding solvers
  KSP* ksp;

  Vec sub_left, coarse_sol;
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
  PetscCall(VecDestroy(&data->coarse_sol));
  PetscCall(PetscFree2(data->ksp, data->mat));
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

PetscErrorCode PCSetup_Net2AS_SetupKSP(PC pc, MPI_Comm comm, PetscInt i) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PC subpc;
  const char* prefix;

  PetscFunctionBegin;
  PetscCall(KSPCreate(comm, data->ksp+i));
  PetscCall(KSPSetType(data->ksp[i], KSPPREONLY));
  PetscCall(KSPGetPC(data->ksp[i], &subpc));
  PetscCall(PCSetType(subpc, PCCHOLESKY));

  PetscCall(PCGetOptionsPrefix(pc, &prefix));
  PetscCall(KSPSetOptionsPrefix(data->ksp[i], prefix));
  PetscCall(KSPAppendOptionsPrefix(data->ksp[i], "net2as_"));
  PetscCall(PCSetOptionsPrefix(subpc, prefix));
  PetscCall(PCAppendOptionsPrefix(subpc, "net2as_"));
  PetscCall(PCSetFromOptions(subpc));
  PetscCall(KSPSetFromOptions(data->ksp[i]));

  PetscCall(KSPSetOperators(data->ksp[i], data->mat[i], data->mat[i]));
  PetscCall(KSPSetUp(data->ksp[i]));
  PetscFunctionReturn(0);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec v = data->vertices;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt vstart, vend, size, nnz = 0, n_cols = data->p[0]*data->p[1], n_rows;
  size_t max_cols;
  PetscInt *rows, *cols;
  PetscReal *vals, h[2], eps = 1e-14;
  PetscBool done;
  const PetscInt *ioff, *inds;
  std::span<PetscReal> vspan;
  Mat coarse_basis, A;
  MatType type;

  PetscFunctionBegin;
  PetscCall(PCDestroy_Net2AS(pc));
  PetscCall(PCGetOperators(pc, &A, NULL));
  PetscCall(VecGetSize(v, &size));
  PetscCall(MatGetSize(A, &data->bs, NULL));
  data->bs /= size/3;
  data->sz = n_cols+1;
  for (PetscInt i = 0; i < 3; i++) {
    PetscCall(VecStrideMin(v, i, NULL, data->min+i));
    PetscCall(VecStrideMax(v, i, NULL, data->max+i));
    if (i < 2) h[i] = (data->max[i]-data->min[i])/(data->p[i]+1);
  }

  PetscCall(VecGetOwnershipRange(v, &vstart, &vend));
  vend /= 3;
  vstart /= 3;
  size = vend-vstart;
  max_cols = data->bs * 4zu * size;
  PetscCall(PetscMalloc3(max_cols, &rows, max_cols, &cols, max_cols, &vals));
  PetscCall(VecGetSpan(v, vspan));
  for (PetscInt n = 0; n < size; n++) {
    PetscReal x = vspan[3*n],            y = vspan[3*n+1];
    PetscInt  i = (x-data->min[0])/h[0], j = (y-data->min[1])/h[1];
    // map to reference element
    PetscReal xx = (x-(i*h[0]+data->min[0]))/h[0], yy = (y-(j*h[1]+data->min[1]))/h[1];

    struct {
      PetscInt row;
      PetscInt col;
      PetscReal val;
    } rcv[4];
    PetscInt n_rcv = 0;

    if (i>data->p[0] || j>data->p[1]) continue;

    if (i>0 && j>0)
      rcv[n_rcv++] = {vstart+n, (j-1)*data->p[0]+(i-1), (1-xx)*(1-yy)};

    if (i<data->p[0] && j>0)
      rcv[n_rcv++] = {vstart+n, (j-1)*data->p[0]+i, xx*(1-yy)};

    if (i>0 && j<data->p[1])
      rcv[n_rcv++] = {vstart+n, j*data->p[0]+i-1, (1-xx)*yy};

    if (i<data->p[0] && j<data->p[1])
      rcv[n_rcv++] = {vstart+n, j*data->p[0]+i, xx*yy};

    for (PetscInt i_rcv = 0; i_rcv < n_rcv; i_rcv++) {
      for (PetscInt i_dofs = 0; i_dofs < data->bs; i_dofs++) {
        rows[nnz] = data->bs * rcv[i_rcv].row + i_dofs;
        cols[nnz] = data->bs * rcv[i_rcv].col + i_dofs;
        vals[nnz] = rcv[i_rcv].val;
        nnz++;
      }
    }
  }
  PetscCall(VecRestoreSpan(v, vspan));

  // setup coarse basis
  PetscCall(MatGetType(A, &type));
  PetscCall(MatCreate(comm, &coarse_basis));
  PetscCall(MatSetType(coarse_basis, type));
  PetscCall(MatSetSizes(coarse_basis, size*data->bs, PETSC_DECIDE, PETSC_DETERMINE, n_cols*data->bs));
  PetscCall(MatSetOptionsPrefix(coarse_basis, "coarse_"));
  PetscCall(MatSetPreallocationCOO(coarse_basis, nnz, rows, cols));
  PetscCall(MatSetValuesCOO(coarse_basis, vals, INSERT_VALUES));
  PetscCall(MatFilter(coarse_basis, eps, /* compress = */ PETSC_TRUE, /* keep = */ PETSC_FALSE));
  PetscCall(PetscFree3(rows, cols, vals));

  // setup sub domains
  PetscCall(MatTranspose(coarse_basis, MAT_INITIAL_MATRIX, &data->sub)); // INITIAL -> redistributes
  PetscCall(MatCreateVecs(data->sub, NULL, &data->sub_left));
  PetscCall(MatCreateVecs(data->sub, NULL, &data->coarse_sol));
  PetscCall(MatViewFromOptions(data->sub, NULL, "-pc_net2as_sub_view"));
  PetscCall(MatGetRowIJ(data->sub, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  PetscCall(MatGetOwnershipRange(data->sub, &vstart, &vend));

  // setup coarse mat
  PetscCall(PetscMalloc2(data->sz, &data->ksp, data->sz, &data->mat));
  PetscCall(MatPtAP(A, coarse_basis, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->mat[0]));
  PetscCall(MatDestroy(&coarse_basis));
  PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_WORLD, 0));

  // setup fine mat
  for (PetscInt s = 0; s < vend-vstart; s++) {
    IS is;
    Mat *mat;
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_USE_POINTER, &is));
    // NOTE: MatCreateSubMatrix creates a submatrix of same type as A, regardless of comm of is,
    //       while MatCreateSubmatrices always creates sequential matrices,
    //       tough it also allocates the output parameter
    PetscCall(MatCreateSubMatrices(A, 1, &is, &is, MAT_INITIAL_MATRIX, &mat));
    data->mat[vstart+s+1] = *mat; // coarse and fine in same array
    PetscCall(PetscFree(mat));
    PetscCall(ISDestroy(&is));
    PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_SELF, vstart+s+1));
  }
  PetscCall(MatRestoreRowIJ(data->sub, 0, PETSC_FALSE, PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  const PetscInt *ioff, *inds;
  PetscInt n_rows, start, end;
  PetscBool done;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));

  // coarse
  PetscCall(MatMult(data->sub, x, data->sub_left)); // sub = coarse_basis^T
  PetscCall(KSPSolve(data->ksp[0], data->sub_left, data->coarse_sol)); // override left // BUG
  PetscCall(MatMultTranspose(data->sub, data->coarse_sol, y));

  PetscCall(MatGetRowIJ(data->sub, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  PetscCall(MatGetOwnershipRange(data->sub, &start, &end));

  // fine
  // TODO: iterate only over the owned rows in sub, do coarse system separately
  for (PetscInt s = start; s < end; s++) {
    PetscInt i = s+1; // subdomain
    const PetscReal *vals;
    IS is;
    Vec z, res;
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_USE_POINTER, &is));
    PetscCall(VecGetSubVector(x, is, &z));
    PetscCall(VecDuplicate(z, &res));
    PetscCall(KSPSolve(data->ksp[i], z, res));
    PetscCall(VecGetArrayRead(res, &vals));
    PetscCall(VecSetValues(y, ioff[s+1]-ioff[s], inds+ioff[s], vals, ADD_VALUES));
    PetscCall(VecRestoreArrayRead(res, &vals));
    PetscCall(VecDestroy(&res));
    PetscCall(VecRestoreSubVector(x, is, &z));
    PetscCall(ISDestroy(&is));
  }

  PetscCall(MatRestoreRowIJ(data->sub, 0, PETSC_FALSE, PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatRestoreRowIJ not done");

  PetscCall(VecAssemblyBegin(y));
  PetscCall(VecAssemblyEnd(y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
  if (!isascii) goto end;

  PetscCall(PetscViewerASCIIPrintf(viewer, "bs=%d\n", data->bs));
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
  PCType type;

  PetscFunctionBeginUser;
  PCGetType(pc, &type);
  if (strcmp(type, "net2as") != 0) PetscFunctionReturn(0);

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
  data->p[0] = data->p[1] = 1;

  pc->ops->apply = PCApply_Net2AS;
  pc->ops->setup = PCSetup_Net2AS;
  pc->ops->destroy = PCDestroy_Net2AS;
  pc->ops->setfromoptions = PCSetFromOptions_Net2AS;
  pc->ops->view = PCView_Net2AS;

  PetscFunctionReturn(PETSC_SUCCESS);
}

template<typename HDG>
PetscErrorCode PCNet2ASVisCoarse(PC pc, HDG& hdg, const char* name) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec left, right;
  PetscInt start, end;
  std::span<PetscReal> span;
  PCType type;

  PetscFunctionBeginUser;
  PCGetType(pc, &type);
  if (strcmp(type, "net2as") != 0) PetscFunctionReturn(0);

  hdg.plot_option("fileName", name);
  hdg.plot_option("printFileNumber", "true");
  PetscCall(MatCreateVecs(data->sub, &right, &left));
  PetscCall(VecGetOwnershipRange(left, &start, &end));
  for (PetscInt i = 0; i < end-start; i++) {
    PetscCall(VecSetValue(left, i+start, 1., INSERT_VALUES));
    PetscCall(MatMultTranspose(data->sub, left, right));
    PetscCall(VecGetSpan(right, span));
    hdg.plot_solution(span, i);
    PetscCall(VecRestoreSpan(right, span));
    PetscCall(VecZeroEntries(left));
  }
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  // constexpr unsigned int poly_deg = 5;
    using Top = Topology::File<1,3>;
    using Geo = Geometry::File<1,3>;
    using NDes = NodeDescriptor::File<1,3>;
    // using LSol = LocalSolver::TimoshenkoBeam<1,3,poly_deg,2*poly_deg,LocalSolver::TimoschenkoBeamParametersClamped>;
    using LSol = LocalSolver::Diffusion<1,5,10,ConstantDiffusionParameters>;
    using HDG = GlobalLoop::Elliptic<Top,Geo,NDes,LSol>;
    constexpr PetscInt bs = LSol::n_glob_dofs_per_node();

    int rank;
    PetscReal rtol = 1e-10;

    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "network";
    char domain_filepath[PATH_MAX] = "domains/grid_8.geo.bin";
    char plot_scale[PATH_MAX] = "1";
    char viscoarse[PATH_MAX] = {0};
    char t2f_mat_load_path[PATH_MAX] = {0};

    PetscLogStage s_as, s_it, s_rf, s_ksp;

    PetscBool is_set, help;
    PetscInt N, ncoo;
    PetscReal tau = 1;
    PetscInt iterations;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    VecScatter scatter;
    Vec rhs, sol, sol0;
    Mat mat;
    KSP ksp;
    PC pc;
    KSPConvergedReason reason;
    PetscReal rnorm;
    PetscBool mat_only = PETSC_FALSE;
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Network Options", NULL);
    PetscCall(PetscOptionsString("-domain", "input network domain", NULL, domain_filepath, domain_filepath, PATH_MAX, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-viscoarse", "output name for visualization of coarse system", NULL, viscoarse, viscoarse, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-mat_only", "only assemble matrix", NULL, mat_only, &mat_only, &is_set));
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

    PetscCall(PetscLogStageRegister("assembly", &s_as));
    PetscCall(PetscLogStageRegister("iteration", &s_it));
    PetscCall(PetscLogStageRegister("residual", &s_rf));
    PetscCall(PetscLogStageRegister("ksp", &s_ksp));

    HDG hdg(domain_filepath);
    zero_v = hdg.zero_vector();
    N = zero_v.size();
    PetscCall(VecCreate(PETSC_COMM_WORLD, &sol));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &rhs));
    PetscCall(VecSetType(sol, VECMPI));
    PetscCall(VecSetType(rhs, VECMPI));
    PetscCall(VecSetSizes(sol, PETSC_DECIDE, N));
    PetscCall(VecSetSizes(rhs, PETSC_DECIDE, N));

    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
    PetscCall(KSPMonitorRegister("csv", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorCSV, NULL, NULL));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(PCNet2ASReadGraph(pc, domain_filepath));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_csv", "csv", NULL));
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", bs, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-t2f_mat_load", t2f_mat_load_path, PATH_MAX, &is_set));
    PetscCall(KSPSetOperators(ksp, mat, mat));

    PRIN2S(s_as);
    if (*t2f_mat_load_path) {
      PetscViewer v;
      PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, t2f_mat_load_path, FILE_MODE_READ, &v));
      PetscCall(MatLoad(mat, v));
    } else {
      mat_coo = hdg.trace_to_flux_mat();
      ncoo = mat_coo.value_vec.size();
      PetscCall(MatSetPreallocationCOO(mat, ncoo, (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
      PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
      PetscCall(MatEliminateZeros(mat, /* keep = */ PETSC_FALSE));
    }
    PRIN2SP();

    if (mat_only) goto end;

    PRIN2S(s_ksp);
    PetscCall(KSPSetUp(ksp));
    PRIN2SP();

    if (strlen(viscoarse) > 0) PetscCall(PCNet2ASVisCoarse(pc, hdg, viscoarse));

    PRIN2S(s_rf);
    PetscCall(VecGetSpan(rhs, rhs_span));
    hdg.residual_flux2(zero_v, rhs_span, 0.);
    PetscCall(VecRestoreSpan(rhs, rhs_span));
    PRIN2SP();

    PetscCall(VecScale(rhs, -1.));

    PRIN2S(s_it);
    PetscCall(KSPSolve(ksp, rhs, sol));
    PRIN2SP();

    PetscCall(KSPGetIterationNumber(ksp, &iterations));
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    PetscCall(KSPSetErrorIfNotConverged(ksp, PETSC_TRUE));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));
    PRIN2IY(iterations);
    PRIN2FY(rnorm);

    PetscCall(VecScatterCreateToZero(sol, &scatter, &sol0));
    PetscCall(VecScatterBegin(scatter, sol, sol0, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, sol, sol0, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecGetSpan(sol0, sol_span));

    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "false");
    hdg.plot_option("scale", plot_scale);
    if (rank == 0) hdg.plot_solution(sol_span);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "output: %s/%s.vtu\n", output_directory, output_filename));

    PetscCall(VecRestoreSpan(sol, sol_span));

end:
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&rhs));
    PetscCall(VecDestroy(&sol));
    PetscCall(VecDestroy(&sol0));
    PetscCall(VecScatterDestroy(&scatter));

    PetscCall(PetscFinalize());
    return 0;
}
