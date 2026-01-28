#include "net2as.hxx"
#include "prin2.hxx"
#include <petsc/private/pcimpl.h>
#include <petscviewerhdf5.h>

struct MatCOO {
  PetscInt *rows, *cols, nnz, cap;
  PetscReal *vals;
};

PetscErrorCode MatCOO_Alloc(MatCOO *coo, PetscInt cap) {
  PetscFunctionBegin;
  PetscCall(PetscMalloc3(cap, &coo->rows, cap, &coo->cols, cap, &coo->vals));
  coo->cap = cap;
  coo->nnz = 0;
  PetscFunctionReturn(0);
};

PetscErrorCode MatCOO_Free(MatCOO *coo) {
  PetscFunctionBegin;
  PetscCall(PetscFree3(coo->rows, coo->cols, coo->vals));
  PetscFunctionReturn(0);
}

PetscErrorCode MatCOO_Push(MatCOO *coo, PetscInt row, PetscInt col, PetscReal val) {
  PetscInt nnz = coo->nnz;

  PetscFunctionBegin;
  PetscAssert(nnz < coo->cap, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "MatCOO_Push out of range");
  coo->rows[nnz] = row;
  coo->cols[nnz] = col;
  coo->vals[nnz] = val;
  coo->nnz++;
  PetscFunctionReturn(0);
}

PetscErrorCode MatCOO_Push_s(MatCOO *coo, PetscInt row, PetscInt col, PetscReal val) {
  PetscInt nnz = coo->nnz;

  PetscFunctionBegin;
  PetscCheck(nnz < coo->cap, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "MatCOO_Push out of range");
  coo->rows[nnz] = row;
  coo->cols[nnz] = col;
  coo->vals[nnz] = val;
  coo->nnz++;
  PetscFunctionReturn(0);
}

struct PC_Net2AS {
  // number of subdomains in [x,y]
  PetscInt p[2];
  // number of local data structures (should be prod(p)+1)
  PetscInt sz;
  // block size (number of dofs per node)
  PetscInt bs;
  // option to print local information
  PetscBool print_local;
  // coarse basis matrix small entry filter tolerance
  PetscReal eps;

  // path to domain file
  char domain[PATH_MAX];
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec points;
  // types 1 -> dirichlet
  IS types_points;
  // dirichlet points
  IS dirichlet;
  // sparse adj matrix representation of the edges in the network
  Mat adj;
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

  PetscCall(PetscOptionsBoundedInt("-net2as_p", "number of subdomains per axis", NULL, p, &p, &set, p_lb));
  if (set) data->p[0] = data->p[1] = p;
  PetscCall(PetscOptionsBoundedInt("-net2as_px", "number of subdomains", NULL, data->p[0], &data->p[0], &set, p_lb));
  PetscCall(PetscOptionsBoundedInt("-net2as_py", "number of subdomains", NULL, data->p[1], &data->p[1], &set, p_lb));
  PetscCall(PetscOptionsBool("-net2as_print_local", "wether to print the local sizes", NULL, data->print_local, &data->print_local, &set));
  PetscCall(PetscOptionsReal("-net2as_eps", "filter tolerance", NULL, data->eps, &data->eps, &set));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscViewer viewer;
  PetscInt m, mm, n, nn, bs;
  PetscInt is_size, is_local, dsize = 0, start, end;
  PetscInt* dir = NULL;
  const PetscInt *types, *ledges;
  IS edges;

  PetscFunctionBeginUser;
  PetscCall(PetscViewerHDF5Open(comm, data->domain, FILE_MODE_READ, &viewer));
  PetscCall(PetscViewerHDF5PushGroup(viewer, "/domain"));
  PetscCall(VecCreate(comm, &data->points));
  PetscCall(PetscObjectSetName((PetscObject)data->points, "points"));
  PetscCall(VecLoad(data->points, viewer));
  PetscCall(VecGetSize(data->points, &n));
  PetscCall(VecGetOwnershipRange(data->points, &start, &end));
  PetscCall(VecGetBlockSize(data->points, &bs));
  end /= bs;
  start /= bs;
  n /= bs;
  nn = end-start;

  PetscCall(ISCreate(comm, &data->types_points));
  PetscCall(PetscObjectSetName((PetscObject)data->types_points, "types_points"));
  PetscCall(ISLoad(data->types_points, viewer));
  PetscCall(ISGetSize(data->types_points, &is_size));
  PetscCall(ISGetLocalSize(data->types_points, &is_local));

  PetscCheck(is_size == n, comm, PETSC_ERR_ARG_SIZ, "types_points size %d != points size %d", is_size, nn);
  PetscCheck(is_local == nn, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "local sizes: types_points %d != points %d", is_local, nn);

  // is_size is local
  PetscCall(PetscMalloc1(is_local, &dir));
  PetscCall(ISGetIndices(data->types_points, &types));
  for (PetscInt i = 0; i < is_local; i++) {
    if (types[i] != 0) dir[dsize++] = start+i;
  }
  PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, dsize, dir, PETSC_OWN_POINTER, &data->dirichlet));
  PetscCall(ISRestoreIndices(data->types_points, &types));

  PetscCall(ISCreate(comm, &edges));
  PetscCall(PetscObjectSetName((PetscObject)edges, "edges"));
  PetscCall(ISLoad(edges, viewer));
  PetscCall(ISGetSize(edges, &m));
  PetscCall(ISGetLocalSize(edges, &mm));
  m /= 2;
  mm /= 2;

  MatCOO coo;
  PetscCall(MatCOO_Alloc(&coo, mm));
  PetscCall(ISGetIndices(edges, &ledges));
  for (PetscInt i = 0; i < mm; i++)
    PetscCall(MatCOO_Push(&coo, ledges[2*i], ledges[2*i+1], 1.));
  PetscCall(ISRestoreIndices(edges, &ledges));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &data->adj));
  PetscCall(MatSetType(data->adj, MATMPIAIJ));
  PetscCall(MatSetSizes(data->adj, m, m, PETSC_DETERMINE, PETSC_DETERMINE));
  PetscCall(MatSetOptionsPrefix(data->adj, "adj_"));
  PetscCall(MatSetPreallocationCOO(data->adj, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(data->adj, coo.vals, INSERT_VALUES));
  PetscCall(MatCOO_Free(&coo));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_domain:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  points: %" PetscInt_FMT "\n", n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  edges: %" PetscInt_FMT "\n", m));
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

PetscErrorCode net2as_cb_q1(PC_Net2AS *data, MatCOO *coo) {
  PetscReal h[2], min[2], max[2];
  PetscInt vstart, vend, size;
  std::span<PetscReal> vspan;

  PetscFunctionBegin;
  for (PetscInt i = 0; i < 2; i++) {
    PetscCall(VecStrideMin(data->points, i, NULL, min+i));
    PetscCall(VecStrideMax(data->points, i, NULL, max+i));
    h[i] = (max[i]-min[i])/(data->p[i]+1);
  }

  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  vend /= 3;
  vstart /= 3;
  size = vend-vstart;
  PetscCall(MatCOO_Alloc(coo, 4zu * size));
  PetscCall(VecGetSpan(data->points, vspan));
  for (PetscInt n = 0; n < size; n++) {
    PetscReal x = vspan[3*n],            y = vspan[3*n+1];
    PetscInt  i = (x-min[0])/h[0], j = (y-min[1])/h[1];
    // map to reference element
    PetscReal xx = (x-(i*h[0]+min[0]))/h[0], yy = (y-(j*h[1]+min[1]))/h[1];

    if (i>data->p[0] || j>data->p[1]) continue;

    if (i>0 && j>0)
      PetscCall(MatCOO_Push(coo, vstart+n, (j-1)*data->p[0]+(i-1), (1-xx)*(1-yy)));
    if (i < data->p[0] && j > 0)
      PetscCall(MatCOO_Push(coo, vstart+n, (j-1)*data->p[0]+i, xx*(1-yy)));
    if (i > 0 && j < data->p[1])
      PetscCall(MatCOO_Push(coo, vstart+n, j*data->p[0]+(i-1), (1-xx)*yy));
    if (i < data->p[0] && j < data->p[1])
      PetscCall(MatCOO_Push(coo, vstart+n, j*data->p[0]+i, xx*yy));
  }
  PetscCall(VecRestoreSpan(data->points, vspan));

  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_alg(PC_Net2AS *data, MatCOO *coo) {
  MatPartitioning p_ctx;
  IS partition, subdoms;
  PetscInt p = data->p[0]*data->p[1];

  PetscFunctionBegin;
  PetscCall(MatPartitioningCreate(PETSC_COMM_WORLD, &p_ctx));
  PetscCall(MatPartitioningSetAdjacency(p_ctx, data->adj));
  PetscCall(MatPartitioningSetNParts(p_ctx, p));
  PetscCall(MatPartitioningSetFromOptions(p_ctx));
  PetscCall(MatPartitioningApply(p_ctx, &partition));
  PetscCall(MatPartitioningDestroy(&p_ctx));

  (void)subdoms;

  // allocate counts
  // count locally
  // create subdom to rank mapping -> this should depend only on the (globally known) number of total subdomains
  // we could assume that the number of subdomains evenly divides the number of procs and assume their sizes to be distributed enough
  // count how many subdoms per rank (maximum)

  PetscFunctionReturn(0);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec gtemp;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt vstart, vend, size, msize, n_cols = data->p[0] * data->p[1], n_rows, n, m;
  PetscBool done;
  const PetscInt *ioff, *inds;
  Mat coarse_basis, A, subdomains;
  MatType type;
  int comm_size;
  MatCOO coo;
  PetscLogDouble t0, t1;

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
  size /= 3;

  PetscCall(net2as_cb_q1(data, &coo));

  // setup coarse basis
  PetscCall(MatGetType(A, &type));
  PetscCall(MatCreate(comm, &coarse_basis));
  PetscCall(MatSetType(coarse_basis, type));
  PetscCall(MatSetSizes(coarse_basis, PETSC_DECIDE, PETSC_DECIDE, size, n_cols));
  PetscCall(MatSetOptionsPrefix(coarse_basis, "coarse_"));
  PetscCall(MatSetPreallocationCOO(coarse_basis, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(coarse_basis, coo.vals, INSERT_VALUES));
  PetscCall(MatZeroRowsIS(coarse_basis, data->dirichlet, 0, NULL, NULL));
  PetscCall(MatFilter(coarse_basis, data->eps, /* compress = */ PETSC_TRUE, /* keep = */ PETSC_FALSE));
  PetscCall(MatCOO_Free(&coo));
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
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  bs: %" PetscInt_FMT "\n", data->bs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  p: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", data->p[0], data->p[1]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  sz: %" PetscInt_FMT "\n", n_cols));
  PetscCall(MatGetSize(A, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  global_size: %" PetscInt_FMT "\n", m));
  PetscCall(MatGetSize(data->mat[0], &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  coarse:\n    size: %" PetscInt_FMT "\n", m));
  PetscCall(PetscTime(&t0));
  PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_WORLD, 0));
  PetscCall(PetscTime(&t1));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    time: %.5e\n", t1-t0));

  if (data->print_local) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  local:\n"));
  PetscCall(MatGetRowIJ(subdomains, 0, /* symmetric = */ PETSC_FALSE, /* inodecomp = */ PETSC_TRUE, &n_rows, &ioff, &inds, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ not done");
  // setup local mat
  for (PetscInt s = 0; s < vend-vstart; s++) {
    Mat *mat;
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, ioff[s+1]-ioff[s], inds+ioff[s], PETSC_COPY_VALUES, &data->is[s+1]));
    if (data->print_local) PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "    - size: %" PetscInt_FMT "\n", data->bs*(ioff[s+1]-ioff[s])));
    // NOTE: MatCreateSubMatrix creates a submatrix of same type as A, regardless of comm of is,
    //       while MatCreateSubmatrices always creates sequential matrices,
    //       tough it also allocates the output parameter
    PetscCall(MatCreateSubMatrices(A, 1, &data->is[s+1], &data->is[s+1], MAT_INITIAL_MATRIX, &mat));
    data->mat[s+1] = *mat;
    PetscCall(PetscFree(mat));
    PetscCall(PetscTime(&t0));
    PetscCall(PCSetup_Net2AS_SetupKSP(pc, PETSC_COMM_SELF, s+1));
    PetscCall(PetscTime(&t1));
    if (data->print_local) PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      time: %.5e\n", t1-t0));
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
  data->eps = 1e-16;

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
