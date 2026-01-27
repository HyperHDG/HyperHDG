#include "net2as.hxx"
#include "prin2.hxx"
#include <petsc/private/pcimpl.h>
#include <petscviewerhdf5.h>

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
  PetscBool print_local_size;

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
  PetscCall(PetscOptionsBool("-pc_net2as_print_local_size", "wether to print the local sizes", NULL, data->print_local_size, &data->print_local_size, &set));
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
  PetscCall(MatGetSize(A, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  global_size: %" PetscInt_FMT "\n", m));
  PetscCall(MatGetSize(data->mat[0], &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  coarse_size: %" PetscInt_FMT "\n", m));
  if (data->print_local_size) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  local_size:\n"));

  PetscCall(MatGetRowIJ(subdomains, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_TRUE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  // setup local mat
  for (PetscInt s = 0; s < vend-vstart; s++) {
    Mat *mat;
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_COPY_VALUES, &data->is[s+1]));
    if (data->print_local_size) PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "    - %" PetscInt_FMT "\n", data->bs*(ioff[s+1]-ioff[s])));
    // NOTE: MatCreateSubMatrix creates a submatrix of same type as A, regardless of comm of is,
    //       while MatCreateSubmatrices always creates sequential matrices,
    //       tough it also allocates the output parameter
    PetscCall(MatCreateSubMatrices(A, 1, &data->is[s+1], &data->is[s+1], MAT_INITIAL_MATRIX, &mat));
    data->mat[s+1] = *mat;
    PetscCall(PetscFree(mat));
    PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_SELF, s+1));
    PetscCall(VecScatterCreate(gtemp, data->is[s+1], data->sol[s+1], NULL, &data->sc[s+1]));
  }
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
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

PetscErrorCode PCNet2ASGetCB(PC pc, Mat *cb) {
  PC_Net2AS *data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  PetscAssertPointer(cb, 2);
  data = (PC_Net2AS *)pc->data;
  *cb = data->cb;
  PetscFunctionReturn(PETSC_SUCCESS);
}
