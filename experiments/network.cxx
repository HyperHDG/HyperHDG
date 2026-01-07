#include <stdio.h>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <petsc/private/pcimpl.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include "parameters.hxx"

static const char help_msg[] = "experiments regarding timoshenko networks\n";
static PetscInt PETSC_PRIN2_ROW_LEN = 10;
static PetscLogDouble PETSC_PRIN2_TIMER = 0;
static PetscInt PETSC_PRIN2_STAGE = 0;
static const char* PETSC_PRIN2_STAGE_NAME = "";

#define PRIN2IY(VAR)  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2FY(VAR)  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2SY(VAR)  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", #VAR, VAR))
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

PetscErrorCode KSPMonitorYAML(KSP ksp, PetscInt it, PetscReal rnorm, PetscViewerAndFormat *vf) {
  PetscViewer viewer = vf->viewer;

  PetscFunctionBegin;
  if (it == 0) PetscCall(PetscViewerASCIIPrintf(viewer, "ksp_monitor:\n"));
  PetscCall(PetscViewerASCIIPrintf(viewer, "  - it: %3" PetscInt_FMT "\n    rnorm: %.16e\n", it, (double)rnorm));
  PetscFunctionReturn(0);
}

struct PC_Net2AS {
  // path to domain file
  char domain[PATH_MAX];
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec points;
  // bounding box of points
  PetscReal min[3], max[3];
  // number of subdomains in [x,y]
  PetscInt p[2];
  // number of local data structures
  PetscInt sz;
  // block size (number of dofs per node)
  PetscInt bs;

  // coarse basis representation of the overlapping subdomains,
  // expanded by block size
  Mat  cb;

  // local datastructures
  // layout: i=0 -> coarse, 1 <= i < sz -> local, corresponding to subdomains
  //   ignore is[0]
  Mat* mat;
  KSP* ksp;
  IS* is;
  Vec* sol;
  VecScatter* sc;
};

PetscErrorCode PCDestroy_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscFunctionBegin;
  for (PetscInt i = 0; data->ksp && i < data->sz; i++) {
    PetscCall(KSPDestroy(data->ksp+i));
    PetscCall(MatDestroy(data->mat+i));
    PetscCall(VecDestroy(data->sol+i));
    if (i > 0) {
      PetscCall(ISDestroy(data->is+i));
      PetscCall(VecScatterDestroy(data->sc+i));
    }
  }
  PetscCall(MatDestroy(&data->cb));
  PetscCall(PetscFree5(data->ksp, data->mat, data->is, data->sol, data->sc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool set;
  PetscInt p = data->p[0], p_lb = data->p[0];

  PetscFunctionBegin;
  PetscCall(PetscOptionsGetString(NULL, NULL, "-domain", data->domain, PATH_MAX, &set));
  PetscOptionsHeadBegin(PetscOptionsObject, "Net2AS options");

  PetscCall(PetscOptionsBoundedInt("-pc_net2as_p", "number of subdomains per axis", NULL, p, &p, &set, p_lb));
  if (set) data->p[0] = data->p[1] = p;
  PetscCall(PetscOptionsBoundedInt("-pc_net2as_px", "number of subdomains", NULL, data->p[0], &data->p[0], &set, p_lb));
  PetscCall(PetscOptionsBoundedInt("-pc_net2as_py", "number of subdomains", NULL, data->p[1], &data->p[1], &set, p_lb));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscViewer viewer;
  PetscInt n, bs;

  PetscFunctionBeginUser;
  PetscCall(PetscViewerHDF5Open(comm, data->domain, FILE_MODE_READ, &viewer));
  PetscCall(PetscViewerHDF5PushGroup(viewer, "/domain"));
  PetscCall(VecCreate(comm, &data->points));
  PetscCall(PetscObjectSetName((PetscObject)data->points, "points"));
  PetscCall(VecLoad(data->points, viewer));
  PetscCall(VecGetSize(data->points, &n));
  PetscCall(VecGetBlockSize(data->points, &bs));
  n /= bs;

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_domain:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, " points: %" PetscInt_FMT "\n", n));
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, " edges: %zu\n", edges.size()));
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

  PetscCall(MatCreateVecs(data->mat[i], &data->sol[i], NULL));
  PetscCall(KSPSetOperators(data->ksp[i], data->mat[i], data->mat[i]));
  PetscCall(KSPSetUp(data->ksp[i]));
  PetscFunctionReturn(0);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec gtemp;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt vstart, vend, size, msize, nnz = 0, n_cols = data->p[0] * data->p[1], n_rows, n, m;
  size_t max_cols;
  PetscInt *rows, *cols;
  PetscReal *vals, h[2], eps = 1e-14;
  PetscBool done;
  const PetscInt *ioff, *inds;
  std::span<PetscReal> vspan;
  Mat coarse_basis, A, subdomains;
  MatType type;
  int comm_size;

  PetscFunctionBegin;

  PetscCallMPI(MPI_Comm_size(comm, &comm_size));
  PetscCheck(n_cols % comm_size == 0, comm, PETSC_ERR_ARG_OUTOFRANGE, "n_cols = %" PetscInt_FMT" must be divisible by MPI_Comm_size = %d", n_cols, comm_size);

  PetscCall(PCDestroy_Net2AS(pc));
  PetscCall(PCSetup_Net2AS_ReadDomain(pc, comm));
  PetscCall(PCGetOperators(pc, &A, NULL));
  PetscCall(VecGetSize(data->points, &size));
  PetscCall(MatGetSize(A, &msize, NULL));
  PetscCall(MatCreateVecs(A, &gtemp, NULL));
  PetscCall(MatGetBlockSize(A, &data->bs));
  PetscCheck(size % 3 == 0, comm, PETSC_ERR_ARG_SIZ,
    "size of points = %" PetscInt_FMT " must be divisible by 3 (x0,y0,z0,x1,y1,z1,...)", size);
  PetscCheck(msize == data->bs * size / 3, comm, PETSC_ERR_ARG_SIZ,
    "A_sz != v_sz / 3 * A_bs where"
    "block size A_bs == %" PetscInt_FMT ", size A_sz == %" PetscInt_FMT ","
    "flat size v_sz == %" PetscInt_FMT, data->bs, msize, size);
  // TODO: take this from hdf5 format
  for (PetscInt i = 0; i < 3; i++) {
    PetscCall(VecStrideMin(data->points, i, NULL, data->min+i));
    PetscCall(VecStrideMax(data->points, i, NULL, data->max+i));
    if (i < 2) h[i] = (data->max[i]-data->min[i])/(data->p[i]+1);
  }

  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  vend /= 3;
  vstart /= 3;
  size = vend-vstart;
  max_cols = 4zu * size;
  PetscCall(PetscMalloc3(max_cols, &rows, max_cols, &cols, max_cols, &vals));
  PetscCall(VecGetSpan(data->points, vspan));
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
      rows[nnz] = rcv[i_rcv].row;
      cols[nnz] = rcv[i_rcv].col;
      vals[nnz] = rcv[i_rcv].val;
      nnz++;
    }
  }
  PetscCall(VecRestoreSpan(data->points, vspan));

  // setup coarse basis
  PetscCall(MatGetType(A, &type));
  PetscCall(MatCreate(comm, &coarse_basis));
  PetscCall(MatSetType(coarse_basis, type));
  PetscCall(MatSetSizes(coarse_basis, size, PETSC_DECIDE, PETSC_DETERMINE, n_cols));
  PetscCall(MatSetOptionsPrefix(coarse_basis, "coarse_"));
  PetscCall(MatSetPreallocationCOO(coarse_basis, nnz, rows, cols));
  PetscCall(MatSetValuesCOO(coarse_basis, vals, INSERT_VALUES));
  PetscCall(MatFilter(coarse_basis, eps, /* compress = */ PETSC_TRUE, /* keep = */ PETSC_FALSE));
  PetscCall(PetscFree3(rows, cols, vals));
  PetscCall(MatCreateMAIJ(coarse_basis, data->bs, &data->cb)); // expanded by block size

  // setup subdomains
  PetscCall(MatTranspose(coarse_basis, MAT_INITIAL_MATRIX, &subdomains)); // INITIAL -> redistributes
  PetscCall(MatViewFromOptions(subdomains, NULL, "-pc_net2as_sub_view"));

  // allocate local datastructures
  PetscCall(MatGetOwnershipRange(subdomains, &vstart, &vend));
  data->sz = 1+vend-vstart;
  PetscCall(PetscMalloc5(data->sz, &data->ksp, data->sz, &data->mat, data->sz, &data->is, data->sz, &data->sol, data->sz, &data->sc));

  // setup coarse mat
  PetscCall(MatPtAP(A, data->cb, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->mat[0]));
  PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_WORLD, 0));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  bs: %" PetscInt_FMT "\n", data->bs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  p: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", data->p[0], data->p[1]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  sz: %" PetscInt_FMT "\n", n_cols));
  (void)m; (void)n;
  // PetscCall(MatGetSize(A, &m, &n));
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  global_size: %" PetscInt_FMT "\n", m));
  // PetscCall(MatGetSize(data->mat[0], &m, &n));
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  coarse_size: %" PetscInt_FMT "\n", m));
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  local_size:\n"));
  // PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  local_sz: %d\n", vend-vstart));
  // PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

  PetscCall(MatGetRowIJ(subdomains, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_TRUE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  // setup local mat
  for (PetscInt s = 0; s < vend-vstart; s++) {
    Mat *mat;
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_COPY_VALUES, &data->is[s+1]));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    - %" PetscInt_FMT "\n", data->bs*(ioff[s+1]-ioff[s])));
    // NOTE: MatCreateSubMatrix creates a submatrix of same type as A, regardless of comm of is,
    //       while MatCreateSubmatrices always creates sequential matrices,
    //       tough it also allocates the output parameter
    PetscCall(MatCreateSubMatrices(A, 1, &data->is[s+1], &data->is[s+1], MAT_INITIAL_MATRIX, &mat));
    data->mat[s+1] = *mat;
    PetscCall(PetscFree(mat));
    PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_SELF, s+1));
    PetscCall(VecScatterCreate(gtemp, data->is[s+1], data->sol[s+1], NULL, &data->sc[s+1]));
  }
  PetscCall(MatRestoreRowIJ(subdomains, 0, PETSC_FALSE, PETSC_FALSE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");

  PetscCall(MatDestroy(&coarse_basis));
  PetscCall(MatDestroy(&subdomains));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));

  // coarse
  PetscCall(MatMultTranspose(data->cb, x, data->sol[0]));
  PetscCall(KSPSolve(data->ksp[0], data->sol[0], data->sol[0]));
  PetscCall(MatMult(data->cb, data->sol[0], y));

  // local
  for (PetscInt i = 1; i < data->sz; i++) {
    PetscCall(VecScatterBegin(data->sc[i], x, data->sol[i], INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(data->sc[i], x, data->sol[i], INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(KSPSolve(data->ksp[i], data->sol[i], data->sol[i]));
    PetscCall(VecScatterBegin(data->sc[i], data->sol[i], y, ADD_VALUES, SCATTER_REVERSE));
    PetscCall(VecScatterEnd(data->sc[i], data->sol[i], y, ADD_VALUES, SCATTER_REVERSE));
  }

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

// NOTE: can only be called sequentially
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
  PetscCall(MatCreateVecs(data->cb, &right, &left));
  PetscCall(VecGetOwnershipRange(left, &start, &end));
  for (PetscInt i = 0; i < end-start; i++) {
    PetscCall(VecSetValue(left, i+start, 1., INSERT_VALUES));
    PetscCall(MatMultTranspose(data->cb, right, left));
    PetscCall(VecGetSpan(left, span));
    hdg.plot_solution(span, i);
    PetscCall(VecRestoreSpan(left, span));
    PetscCall(VecZeroEntries(right));
  }
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
    constexpr unsigned int poly_deg = 5;
    using Top = Topology::File<1,3>;
    using Geo = Geometry::File<1,3>;
    using NDes = NodeDescriptor::File<1,3>;
    using LSol = LocalSolver::TimoshenkoBeam<1,3,poly_deg,2*poly_deg,LocalSolver::TimoschenkoBeamParametersClamped>;
    // using LSol = LocalSolver::Diffusion<1,5,10,ConstantDiffusionParameters>;
    using HDG = GlobalLoop::Elliptic<Top,Geo,NDes,LSol>;
    constexpr PetscInt bs = LSol::n_glob_dofs_per_node();

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

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# ------- " __FILE__ " -------\n"));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
    PetscCallMPI(MPI_Get_processor_name(proc_name, &proc_name_len));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Network Options", NULL);
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

    HDG hdg(domain_filepath);
    PetscCall(PCRegister("net2as", PCCreate_Net2AS));
    PetscCall(KSPMonitorRegister("yaml", PETSCVIEWERASCII, PETSC_VIEWER_DEFAULT, KSPMonitorYAML, NULL, NULL));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, "net2as"));
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPMonitorSetFromOptions(ksp, "-ksp_monitor_yaml", "yaml", NULL));
    PetscCall(KSPSetFromOptions(ksp));

    N = hdg.size_of_system();
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
      mat_coo = hdg.trace_to_flux_mat();
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
      hdg.residual_flux2(hdg.zero_vector(), span, 0.);
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

    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "false");
    hdg.plot_option("scale", plot_scale);
    if (rank == 0) {
      PetscCall(VecGetSpan(rhs0, span));
      hdg.plot_solution(span);
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
