#include "net2as.hxx"
#include "prin2.hxx"
#include <petsc/private/pcimpl.h>
#include <petsc/private/hashmapi.h>
#include <petsc/private/hashseti.h>
#include <petscviewerhdf5.h>

#ifdef HYPERHDG_PARHIP
  #include "petsc_parhip.h"
#endif

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

PetscErrorCode MatCOO_Realloc(MatCOO *coo, PetscInt cap) {
  MatCOO old = *coo;

  PetscFunctionBegin;
  PetscCall(PetscMalloc3(cap, &coo->rows, cap, &coo->cols, cap, &coo->vals));
  PetscCall(PetscArraycpy(old.rows, coo->rows, coo->nnz));
  PetscCall(PetscArraycpy(old.cols, coo->cols, coo->nnz));
  PetscCall(PetscArraycpy(old.vals, coo->vals, coo->nnz));
  coo->cap = cap;
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

PetscErrorCode MatCOO_View(MatCOO *coo, PetscViewer viewer) {
  int rank;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCall(PetscViewerASCIIPrintf(viewer, "MatCOO\n"));
  PetscCall(PetscViewerASCIIPushSynchronized(viewer));
  PetscCall(PetscViewerASCIISynchronizedPrintf(viewer, "  [%d] nnz=%d\n", rank, coo->nnz));
  for (PetscInt i = 0; i < coo->nnz; i++) {
    PetscCall(PetscViewerASCIISynchronizedPrintf(viewer, "  [%d] %6d %6d %.5e\n", rank, coo->rows[i], coo->cols[i], coo->vals[i]));
  }
  PetscCall(PetscViewerASCIIPopSynchronized(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

struct PC_Net2AS {
  // configuration paramters

  // number of subdomains in [x,y]
  PetscInt p[2];
  // number of global coarse dofs
  PetscInt n_coarse;
  // number of local data structures
  PetscInt sz;
  // block size (number of dofs per node)
  PetscInt bs;
  // option to print local information
  PetscBool print_local;
  // coarse basis matrix small entry filter tolerance
  PetscReal eps;
  // bound on the (pointwise) multiplicity of the cover formed by the subdomains
  PetscInt mult_bound;
  // overlap parameter in number of hops
  PetscInt delta;
  // the number of dimension to extend the PU to
  PetscInt pux_dim;
  // cb_type one of "q1", "pu"
  char cb_type[10];
  // load_type one of "rr", "gr"
  // rr - naive round robin load balancing
  // gr - simplest greedy load balancing
  char load_type[10];

  // network information

  // path to domain file
  char domain[PATH_MAX];
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec points;
  // boundary points
  IS boundary;
  // partition of points (only used by cb_pu), same parallel layout as points and as adj
  IS partition;
  // sparse adj matrix representation of the edges in the network
  Mat adj;

  // coarse global data structures
  Mat cb;
  Mat cmat;
  KSP cksp;
  Vec csol;

  // local rank data structures
  VecScatter rank_sc;
  Vec rank_sol;
  IS rank_is;

  // local subdomain datastructures,
  // arrays of length data->sz
  PetscInt *sd_gids; // local subdomain global ids
  Mat* mat;
  KSP* ksp;
  IS* is;
  IS* local_is;
  Vec* sol;
  VecScatter* sc;

  PetscReal bal_est;
};

PetscErrorCode net2as_alloc_ds(PC_Net2AS *data, PetscInt sz) {
  PetscInt min_sz = PetscMax(1, sz);
  PetscFunctionBegin;
  PetscCall(PetscMalloc6(min_sz, &data->ksp, min_sz, &data->is, min_sz,
    &data->sol, min_sz, &data->sc, min_sz, &data->sd_gids, min_sz, &data->local_is));
  data->sz = min_sz;
  PetscFunctionReturn(0);
}

PetscErrorCode PCDestroy_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscFunctionBegin;

  // global resources
  PetscCall(MatDestroy(&data->cb));
  PetscCall(MatDestroy(&data->cmat));
  PetscCall(KSPDestroy(&data->cksp));
  PetscCall(VecDestroy(&data->csol));

  // local rank resources
  PetscCall(VecScatterDestroy(&data->rank_sc));
  PetscCall(VecDestroy(&data->rank_sol));
  PetscCall(ISDestroy(&data->rank_is));

  // local subdomain resources
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(KSPDestroy(data->ksp+i));
    PetscCall(VecDestroy(data->sol+i));
    PetscCall(ISDestroy(data->is+i));
    PetscCall(ISDestroy(data->local_is+i));
  }
  PetscCall(MatDestroySubMatrices(data->sz, &data->mat));
  PetscCall(PetscFree6(data->ksp, data->is, data->sol, data->sc, data->sd_gids, data->local_is));
  PetscCall(ISDestroy(&data->partition));

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
  PetscCall(PetscOptionsInt("-net2as_mult", "upper bound on the (pointwise) multiplicity of the cover formed by the subdomains", NULL, data->mult_bound, &data->mult_bound, &set));
  PetscCall(PetscOptionsInt("-net2as_delta", "overlap parameters in number of hops", NULL, data->delta, &data->delta, &set));
  PetscCall(PetscOptionsInt("-net2as_pux_dim", "number of dimension to extend the pu by", NULL, data->pux_dim, &data->pux_dim, &set));
  PetscCall(PetscOptionsString("-net2as_cb_type", "subdomain partition type", NULL, data->cb_type, data->cb_type, sizeof(data->cb_type), &set));
  PetscCall(PetscOptionsString("-net2as_load_type", "subdomain load balancing type", NULL, data->load_type, data->load_type, sizeof(data->load_type), &set));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscBool net2as_is_dirichlet(PetscInt type, PetscBool wave) {
  if (type == 0) return PETSC_FALSE;
  if (wave && (type & (1u << 6))) return PETSC_FALSE;  // static-only, free in wave
  return PETSC_TRUE;
}

PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscViewer viewer;
  PetscInt m, mm, n, nn, bs;
  PetscInt start, end;
  const PetscInt *ledges;
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

  PetscCall(ISCreate(comm, &edges));
  PetscCall(PetscObjectSetName((PetscObject)edges, "edges"));
  PetscCall(ISLoad(edges, viewer));
  PetscCall(ISGetSize(edges, &m));
  PetscCall(ISGetLocalSize(edges, &mm));
  m /= 2;
  mm /= 2;

  MatCOO coo;
  PetscCall(MatCOO_Alloc(&coo, 2*mm));
  PetscCall(ISGetIndices(edges, &ledges));
  for (PetscInt i = 0; i < mm; i++) {
    PetscCall(MatCOO_Push(&coo, ledges[2*i], ledges[2*i+1], 1.));
    PetscCall(MatCOO_Push(&coo, ledges[2*i+1], ledges[2*i], 1.));
  }
  PetscCall(ISRestoreIndices(edges, &ledges));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &data->adj));
  PetscCall(MatSetType(data->adj, MATMPIAIJ));
  PetscCall(MatSetSizes(data->adj, nn, nn, n, n));
  PetscCall(MatSetOptionsPrefix(data->adj, "net2as_adj_"));
  PetscCall(MatSetPreallocationCOO(data->adj, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(data->adj, coo.vals, INSERT_VALUES));
  PetscCall(MatCOO_Free(&coo));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_domain:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  filepath: %s\n", data->domain));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  points: %" PetscInt_FMT "\n", n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  edges: %" PetscInt_FMT "\n", m));
  PetscFunctionReturn(PETSC_SUCCESS);
}

struct Net2AS_SolveInfo {
  PetscInt size;
  PetscInt nz_mat;
  PetscInt nz_fac;
  PetscLogDouble fill;
  PetscLogDouble time;
  MatSolverType factor_type;
};

PetscErrorCode net2as_setup_ds(PC pc, MPI_Comm comm, const char *prefix, KSP *ksp, Mat *mat, Vec *sol, Net2AS_SolveInfo *info) {
  PC subpc;
  PetscLogDouble t0, t1;
  Mat factored;
  MatInfo minfo;

  PetscFunctionBegin;
  PetscCall(KSPCreate(comm, ksp));
  PetscCall(KSPSetType(*ksp, KSPPREONLY));
  PetscCall(KSPGetPC(*ksp, &subpc));
  PetscCall(PCSetType(subpc, PCCHOLESKY));

  PetscCall(KSPSetOptionsPrefix(*ksp, "net2as_"));
  PetscCall(KSPAppendOptionsPrefix(*ksp, prefix));
  PetscCall(PCSetOptionsPrefix(subpc, "net2as_"));
  PetscCall(PCAppendOptionsPrefix(subpc, prefix));
  PetscCall(PCSetFromOptions(subpc));
  PetscCall(KSPSetFromOptions(*ksp));

  PetscCall(MatCreateVecs(*mat, sol, NULL));
  PetscCall(KSPSetOperators(*ksp, *mat, *mat));
  PetscCall(PetscTime(&t0));
  PetscCall(KSPSetUp(*ksp));
  PetscCall(PetscTime(&t1));
  if (info) {
    info->time = t1-t0;
    PetscCall(MatGetSize(*mat, &info->size, NULL));
    PetscCall(PCFactorGetMatrix(subpc, &factored));
    PetscCall(MatGetInfo(*mat, MAT_GLOBAL_SUM, &minfo));
    info->nz_mat = (PetscInt)minfo.nz_used;
    PetscCall(MatGetInfo(factored, MAT_GLOBAL_SUM, &minfo));
    info->nz_fac = (PetscInt)minfo.nz_used;
    info->fill = (PetscLogDouble)info->nz_fac / info->nz_mat;
    PetscCall(PCFactorGetMatSolverType(subpc, &info->factor_type));
  }

  PetscFunctionReturn(0);
}

// simplest of all greedy load balancing strategies
PetscErrorCode net2as_loadbalance_greedy(MPI_Comm comm, PetscInt *weights, PetscInt *assignments, PetscInt count, PetscReal *bal) {
  int size;
  PetscHeap loads;
  PetscInt max_load = 0, sum_load = 0;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCall(PetscHeapCreate(size, &loads));
  for (PetscInt r = 0; r < size; r++)
    PetscCall(PetscHeapAdd(loads, r, 0));
  for (PetscInt i = 0; i < count; i++) {
    PetscInt r, load;
    PetscCall(PetscHeapPop(loads, &r, &load));
    assignments[i] = r;
    load += weights[i];
    sum_load += weights[i];
    max_load = PetscMax(max_load, load);
    PetscCall(PetscHeapAdd(loads, r, load));
  }
  *bal = (PetscReal)max_load / sum_load * size;
  PetscCall(PetscHeapDestroy(&loads));
  PetscFunctionReturn(0);
}

// naive round robin scheduling
PetscErrorCode net2as_loadbalance_round_robin(MPI_Comm comm, PetscInt *weights, PetscInt *assignments, PetscInt count, PetscReal *bal) {
  int size;
  PetscInt *loads, max_load = 0, sum_load = 0;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCall(PetscMalloc1(size, &loads));
  for (PetscInt i = 0; i < count; i++) {
    PetscInt r = i % size;
    assignments[i] = r;
    loads[r] += weights[i];
    max_load = PetscMax(max_load, weights[i]);
    sum_load += weights[i];
  }
  *bal = (PetscReal)max_load / sum_load * size;
  PetscCall(PetscFree(loads));
  PetscFunctionReturn(0);
}

// cb = coarse basis
//   rows = global ids of local vertices
//   cols = global subdom ids
// sd = subdomains
//   rows = global ids of local subdoms
//   cols = global vertex ids
//   will be allocated, must be freed
PetscErrorCode net2as_distribute_subdomains(MPI_Comm comm, PC_Net2AS *data, MatCOO *cb, MatCOO *sd) {
  int rank, size, tag_vid = 0, tag_sid = 1;
  PetscInt p = data->n_coarse, sd_count, sd_total_size, off, start;
  PetscInt *sd2lcounts, *sd2gcounts, *sd2rank, *rank2scount, *rank2rcount, *coo2rank;
  MPI_Request *reqs;
  PetscInt *tmp_rows, *tmp_cols;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCall(PetscCalloc7(p, &sd2lcounts, p, &sd2gcounts, p, &sd2rank,
    size, &rank2scount, size, &rank2rcount, cb->nnz, &coo2rank, 4*size, &reqs));
  PetscCall(PetscMalloc2(cb->nnz, &tmp_rows, cb->nnz, &tmp_cols));

  // PetscCall(MatCOO_View(cb, PETSC_VIEWER_STDOUT_WORLD));

  // compute sizes of local parts of the subdomains
  for (PetscInt i = 0; i < cb->nnz; i++) sd2lcounts[cb->cols[i]]++;

  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "sd2lcounts", sd2lcounts, p));

  // compute global sizes of the subdomains
  PetscCallMPI(MPI_Allreduce(sd2lcounts, sd2gcounts, p, MPIU_INT, MPI_SUM, comm));

  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "sd2gcounts", sd2gcounts, p));

  // compute some load balancing strategy
  // sd2rank is an assignment of subdomains (indices) to ranks (values)
  if (rank == 0) {
    if (strcmp(data->load_type, "rr") == 0)
      PetscCall(net2as_loadbalance_round_robin(comm, sd2gcounts, sd2rank, p, &data->bal_est));
    else if (strcmp(data->load_type, "gr") == 0)
      PetscCall(net2as_loadbalance_greedy(comm, sd2gcounts, sd2rank, p, &data->bal_est));
    else
      PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
        "unexptected load balancing type '%s', expected one of 'rr', 'gr'", data->load_type);
  }
  // send this assignment to all ranks
  PetscCallMPI(MPI_Bcast(sd2rank, p, MPIU_INT, 0, comm));

  // count how many vertices I will send to each proc
  sd_count = 0;
  sd_total_size = 0;
  for (PetscInt i = 0; i < p; i++) {
    if (sd2rank[i] == rank) {
      sd_count++;
      sd_total_size += sd2gcounts[i];
    }
    rank2scount[sd2rank[i]] += sd2lcounts[i];
  }

  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "sd_count", &sd_count, 1));
  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "sd_total_size", &sd_total_size, 1));
  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "rank2scount", rank2scount, size));

  PetscCall(net2as_alloc_ds(data, sd_count));

  // NOTE: could use petsc build two sided
  // now rank2rcount specifies from which remote rank the local rank will receive how many vertex, sd_id pairs
  PetscCallMPI(MPI_Alltoall(rank2scount, 1, MPIU_INT, rank2rcount, 1, MPIU_INT, comm));

  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "rank2rcount", rank2rcount, size));

  // allocate enough space for the sd the local rank owns
  PetscCall(MatCOO_Alloc(sd, sd_total_size));

  // setup the receives
  // corresponding vertex and subdomain ids
  // NOTE: should first setup receives, otherwise might run into buffering issues
  off = 0;
  for (PetscInt r = 0; r < size; r++) {
    PetscCallMPI(MPI_Irecv(sd->cols+off, rank2rcount[r], MPIU_INT, r, tag_vid, comm, &reqs[2*r]));
    PetscCallMPI(MPI_Irecv(sd->rows+off, rank2rcount[r], MPIU_INT, r, tag_sid, comm, &reqs[2*r+1]));
    off += rank2rcount[r];
  }
  sd->nnz = off;

  // sort cb by ranks of subdomain indices
  for (PetscInt i = 0; i < cb->nnz; i++) {
    coo2rank[i] = sd2rank[cb->cols[i]];
    tmp_rows[i] = cb->rows[i];
    tmp_cols[i] = cb->cols[i];
  }
  PetscCall(PetscSortIntWithArrayPair(cb->nnz, coo2rank, tmp_rows, tmp_cols));

  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "cb sorted by sd2rank\n"));
  // PetscCall(MatCOO_View(cb, PETSC_VIEWER_STDOUT_WORLD));

  // setup sends
  // corresponding vertex and subdomain ids
  off = 0;
  for (PetscInt r = 0; r < size; r++) {
    PetscCallMPI(MPI_Isend(cb->rows+off, rank2scount[r], MPIU_INT, r, tag_vid, comm, &reqs[2*(size+r)]));
    PetscCallMPI(MPI_Isend(cb->cols+off, rank2scount[r], MPIU_INT, r, tag_sid, comm, &reqs[2*(size+r)+1]));
    off += rank2scount[r];
  }

  // then wait on all receives and sends
  PetscCall(MPI_Waitall(4*size, reqs, MPI_STATUSES_IGNORE));

  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "sd unsorted\n"));
  // PetscCall(MatCOO_View(sd, PETSC_VIEWER_STDOUT_WORLD));

  // sort sd by subdomain indices
  PetscCall(PetscSortIntWithArray(sd->nnz, sd->rows, sd->cols));

  // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "sd sorted by subdom\n"));
  // PetscCall(MatCOO_View(sd, PETSC_VIEWER_STDOUT_WORLD));

  // create IS
  start = 0; // start of contiguous subdomain indices
  off = 0;   // local subdomain index in local ds
  while (start < sd->nnz) {
    PetscInt end = start;
    PetscInt *global_vertex_ids = &sd->cols[start];
    while (end < sd->nnz && sd->rows[end] == sd->rows[start])
      end++;
    data->sd_gids[off] = sd->rows[start];
    // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "is", global_vertex_ids, end-start));
    // PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, end-start, global_vertex_ids, PETSC_COPY_VALUES, &data->is[off]));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, end-start, global_vertex_ids, PETSC_COPY_VALUES, &data->is[off]));
    start = end;
    off++;
  }
  PetscCheck(off == sd_count, PETSC_COMM_SELF, PETSC_ERR_PLIB, "detected '%d' subdomains, expected '%d'", off, sd_count);

  if (sd_count == 0 && data->sz > 0)
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, 0, NULL, PETSC_COPY_VALUES, &data->is[0]));

  PetscCall(PetscFree2(tmp_rows, tmp_cols));
  PetscCall(PetscFree7(sd2lcounts, sd2gcounts, sd2rank,
    rank2scount, rank2rcount, coo2rank, reqs));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_single_is_blocked(IS *is, PetscInt bs) {
  const PetscInt *inds;
  PetscInt sz;
  IS bis;

  PetscFunctionBegin;
    PetscCall(ISGetIndices(*is, &inds));
    PetscCall(ISGetLocalSize(*is, &sz));
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, bs, sz, inds, PETSC_COPY_VALUES, &bis));
    PetscCall(ISRestoreIndices(*is, &inds));
    PetscCall(ISDestroy(is));
    *is = bis;

  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_is_blocked(PC_Net2AS *data) {
  PetscFunctionBegin;
  PetscCall(net2as_make_single_is_blocked(&data->rank_is, data->bs));
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(net2as_make_single_is_blocked(data->is+i, data->bs));
    PetscCall(net2as_make_single_is_blocked(data->local_is+i, data->bs));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_rank_is(IS *is, PetscInt sz, PetscInt bs, IS *ris) {
  PetscInt *all, n = 0, nn, off = 0;
  const PetscInt *inds;

  PetscFunctionBegin;
  (void)bs;
  for (PetscInt i = 0; i < sz; i++) {
    PetscCall(ISGetLocalSize(is[i], &nn));
    n += nn;
  }
  PetscCall(PetscMalloc1(n, &all));
  for (PetscInt i = 0; i < sz; i++) {
    PetscCall(ISGetLocalSize(is[i], &nn));
    PetscCall(ISGetIndices(is[i], &inds));
    PetscCall(PetscArraycpy(all+off, inds, nn));
    PetscCall(ISRestoreIndices(is[i], &inds));
    off += nn;
  }
  PetscCall(PetscSortRemoveDupsInt(&off, all));
  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, off, all, PETSC_COPY_VALUES, ris));
  PetscCall(PetscFree(all));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_is_local(PC_Net2AS *data, MatCOO *sd) {
  ISLocalToGlobalMapping l2g;

  PetscFunctionBegin;
  PetscCall(ISLocalToGlobalMappingCreateIS(data->rank_is, &l2g));
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscInt n_global, n_local;

    // NOTE: the new IS is not a block-IS
    PetscCall(ISGlobalToLocalMappingApplyIS(l2g, IS_GTOLM_DROP, data->is[i], &data->local_is[i]));
    PetscCall(ISGetLocalSize(data->is[i], &n_global));
    PetscCall(ISGetLocalSize(data->local_is[i], &n_local));
    PetscCheck(n_global == n_local, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
      "rank_is local to global mapping inds dropped: expected %" PetscInt_FMT ", got %" PetscInt_FMT, n_global, n_local);
  }
  PetscCall(ISLocalToGlobalMappingDestroy(&l2g));

  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_q1(PC_Net2AS *data, MatCOO *coo, MatCOO *sd) {
  PetscReal min[2], max[2], h[2];
  PetscInt vstart, vend, ns[2];
  PetscReal *points;
  PetscInt n_coarse, n_local;
  PetscFunctionBegin;

  for (PetscInt d = 0; d < 2; d++) {
    PetscCall(VecStrideMin(data->points, d, NULL, min+d));
    PetscCall(VecStrideMax(data->points, d, NULL, max+d));
    h[d] = (max[d]-min[d]) / (data->p[d]+1);
    ns[d] = data->p[d]+2; // coarse DoFs per dim
  }
  n_coarse = ns[0] * ns[1];

  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  n_local = (vend - vstart) / 3;
  vstart /= 3;

  PetscCall(VecGetArray(data->points, &points));
  data->n_coarse = n_coarse;

  // fill coarse basis functions
  PetscCall(MatCOO_Alloc(coo, 4 * n_local));
  for (PetscInt k = 0; k < n_local; k++) {
    PetscReal x = points[3*k], y = points[3*k+1];
    PetscInt i = (x - min[0]) / h[0];
    PetscInt j = (y - min[1]) / h[1];
    if (i > data->p[0]) i = data->p[0];
    if (j > data->p[1]) j = data->p[1];
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    PetscReal xx = (x - (i*h[0] + min[0])) / h[0];
    PetscReal yy = (y - (j*h[1] + min[1])) / h[1];

    PetscInt row = vstart + k;
    struct { PetscInt i, j; PetscReal w; } pts[4] = {
      {i,   j,   (1-xx)*(1-yy)},
      {i+1, j,       xx*(1-yy)},
      {i,   j+1, (1-xx)*yy    },
      {i+1, j+1,     xx*yy    },
    };
    for (unsigned int l = 0; l < 4; l++) {
      PetscInt col = pts[l].j * ns[0] + pts[l].i;
      PetscCall(MatCOO_Push(coo, row, col, pts[l].w));
      // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%.2e, %.2e, %d, %d, | %d, %d, %.2e\n", x, y, i, j , row, col, pts[l].w));
    }
  }


  PetscCall(VecRestoreArray(data->points, &points));

  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, coo, sd));
  PetscCall(net2as_make_rank_is(data->is, data->sz, data->bs, &data->rank_is));
  PetscCall(net2as_make_is_local(data, sd));
  PetscCall(net2as_make_is_blocked(data));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_pu(PC_Net2AS *data, MatCOO *coo, MatCOO *sd) {
  MatPartitioning p_ctx;
  MatPartitioningType p_type;
  IS partition = data->partition, bis;
  PetscInt p = data->p[0]*data->p[1], lsz_part, new_cap, vstart, vend, cut, *counts;
  const PetscInt *inds;
  Vec rank_points;
  VecScatter sc;
  PetscReal *points;

  PetscFunctionBegin;

  data->n_coarse = p;
#ifdef HYPERHDG_PARHIP
  PetscCall(MatPartitioningRegister("parhip", MatPartitioningCreate_ParHIP));
#endif

  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  vstart /= 3; vend /= 3;
  PetscCall(MatPartitioningCreate(PETSC_COMM_WORLD, &p_ctx));
  PetscCall(MatPartitioningSetAdjacency(p_ctx, data->adj));
  PetscCall(MatPartitioningSetNParts(p_ctx, p));
  PetscCall(MatPartitioningSetFromOptions(p_ctx));
  PetscCall(MatPartitioningApply(p_ctx, &partition));
  PetscCall(MatPartitioningParmetisGetEdgeCut(p_ctx, &cut));
  PetscCall(MatPartitioningGetType(p_ctx, &p_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_cb_pu:\n  cut: %" PetscInt_FMT "\n", cut));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  pux_dim: %" PetscInt_FMT "\n", data->pux_dim));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  part_type: %s\n", p_type));
  PetscCall(MatPartitioningDestroy(&p_ctx));

  PetscCall(ISGetIndices(partition, &inds));
  PetscCall(ISGetLocalSize(partition, &lsz_part));
  PetscCheck(vend-vstart == lsz_part, PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
    "parallel layout of data->points must match that of the partition, data->points size '%d', "
    "partition size '%d'", vend-vstart, lsz_part);
  PetscCall(MatCOO_Alloc(coo, lsz_part));
  for (PetscInt i = 0; i < lsz_part; i++)
    PetscCall(MatCOO_Push(coo, vstart+i, inds[i], 1.));
  PetscCall(ISRestoreIndices(partition, &inds));

  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, coo, sd));

  PetscCall(MatIncreaseOverlap(data->adj, data->sz, data->is, data->delta));

  PetscCall(net2as_make_rank_is(data->is, data->sz, data->bs, &data->rank_is));
  PetscCall(net2as_make_is_local(data, sd));

  new_cap = 0;
  PetscCall(ISGetLocalSize(data->rank_is, &new_cap));
  PetscCall(PetscCalloc1(new_cap, &counts));
  new_cap = 0;
  for (PetscInt s = 0; s < data->sz; s++) {
    PetscInt sz;
    PetscCall(ISGetIndices(data->local_is[s], &inds));
    PetscCall(ISGetLocalSize(data->local_is[s], &sz));
    for (PetscInt i = 0; i < sz; i++) counts[inds[i]]++;
    PetscCall(ISRestoreIndices(data->local_is[s], &inds));
    new_cap += sz;
  }

  PetscCall(MatCOO_Free(coo));
  PetscCall(MatCOO_Alloc(coo, new_cap*(data->pux_dim+1)));

  PetscCall(ISGetLocalSize(data->rank_is, &new_cap));
  PetscCall(ISGetIndices(data->rank_is, &inds));
  PetscCall(ISCreateBlock(PETSC_COMM_SELF, 3, new_cap, inds, PETSC_COPY_VALUES, &bis));
  PetscCall(ISRestoreIndices(data->rank_is, &inds));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, 3*new_cap, &rank_points));
  PetscCall(VecScatterCreate(data->points, bis, rank_points, NULL, &sc));
  PetscCall(VecScatterBegin(sc, data->points, rank_points, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(sc, data->points, rank_points, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterDestroy(&sc));
  PetscCall(ISDestroy(&bis));
  PetscCall(VecGetArray(rank_points, &points));

  for (PetscInt s = 0; s < data->sz; s++) {
    PetscInt sz, sz2;
    const PetscInt *inds_l;
    PetscReal centroid[3] = {0}, min_sd[3] = {PETSC_MAX_REAL, PETSC_MAX_REAL, PETSC_MAX_REAL}, max_sd[3] = {PETSC_MIN_REAL, PETSC_MIN_REAL, PETSC_MIN_REAL};
    PetscCall(ISGetIndices(data->local_is[s], &inds_l));
    PetscCall(ISGetIndices(data->is[s], &inds));
    PetscCall(ISGetLocalSize(data->is[s], &sz));
    PetscCall(ISGetLocalSize(data->local_is[s], &sz2));
    PetscCheck(sz == sz2, PETSC_COMM_SELF, PETSC_ERR_PLIB, "local and global is size should be the same, %" PetscInt_FMT " != %" PetscInt_FMT, sz, sz2);
    for (PetscInt i = 0; i < sz; i++) {
      for (PetscInt j = 0; j < 3; j++) {
        PetscReal r = points[3*inds_l[i]+j];
        centroid[j] += r;
        min_sd[j] = PetscMin(min_sd[j], r);
        max_sd[j] = PetscMax(max_sd[j], r);
      }
    }
    for (PetscInt j = 0; j < 3; j++) centroid[j] /= sz;

    for (PetscInt i = 0; i < sz; i++) {
      PetscInt li = inds_l[i], cb_idx = (data->pux_dim+1)*data->sd_gids[s];
      PetscCall(MatCOO_Push(coo, inds[i], cb_idx, 1./counts[li]));
      for (PetscInt j = 0; j < data->pux_dim; j++)
        PetscCall(MatCOO_Push(coo, inds[i], cb_idx+j+1, (points[3*li+j]-centroid[j])/(max_sd[j]-min_sd[j]) /counts[li]));
    }
    PetscCall(ISRestoreIndices(data->is[s], &inds));
    PetscCall(ISRestoreIndices(data->local_is[s], &inds_l));
  }

  PetscCall(VecRestoreArray(rank_points, &points));
  PetscCall(VecDestroy(&rank_points));
  PetscCall(PetscFree(counts));

  PetscCall(net2as_make_is_blocked(data));

  PetscFunctionReturn(0);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec gtemp;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt size, msize, n_cols, n, m;
  Mat coarse_basis, A;
  MatType type;
  MatCOO coo, sd;
  Net2AS_SolveInfo info;
  MatInfo mat_info;
  PetscReal t_max = 0, t_sum = 0, t_loc = 0;

  PetscFunctionBegin;

  PetscCall(PCDestroy_Net2AS(pc));
  PetscCallMPI(MPI_Comm_size(comm, &size));
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

  if (strcmp(data->cb_type, "q1") == 0)
    PetscCall(net2as_cb_q1(data, &coo, &sd));
  else if (strcmp(data->cb_type, "pu") == 0)
    PetscCall(net2as_cb_pu(data, &coo, &sd));
  else
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "unsupported type '%s', muse be one of 'q1', 'pu'", data->cb_type);

  n_cols = data->n_coarse;

  // setup coarse basis
  PetscCall(MatGetType(A, &type));
  PetscCall(MatCreate(comm, &coarse_basis));
  PetscCall(MatSetType(coarse_basis, type));
  PetscCall(MatSetSizes(coarse_basis, PETSC_DECIDE, PETSC_DECIDE, size, n_cols*(data->pux_dim+1)));
  PetscCall(MatSetOptionsPrefix(coarse_basis, "net2as_coarse_"));
  PetscCall(MatSetPreallocationCOO(coarse_basis, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(coarse_basis, coo.vals, INSERT_VALUES));
  PetscCall(MatCreateMAIJ(coarse_basis, data->bs, &data->cb)); // expanded by block size

  // setup coarse global data structures
  PetscCall(MatPtAP(A, data->cb, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->cmat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  cb_type: %s\n", data->cb_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  load_type: %s\n", data->load_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  bs: %" PetscInt_FMT "\n", data->bs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  p: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", data->p[0], data->p[1]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  sz: %" PetscInt_FMT "\n", n_cols));
  PetscCall(MatGetSize(A, &m, &n));
  PetscCall(MatGetInfo(A, MAT_GLOBAL_SUM, &mat_info));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  global:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    size: %" PetscInt_FMT "\n", m));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    nz_mat: %" PetscInt_FMT "\n", (PetscInt)mat_info.nz_used));
  PetscCall(MatGetSize(data->cmat, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  coarse:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    size: %" PetscInt_FMT "\n", m));
  PetscCall(net2as_setup_ds(pc, PETSC_COMM_WORLD, "coarse_", &data->cksp, &data->cmat, &data->csol, &info));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    nz_mat: %" PetscInt_FMT "\n", info.nz_mat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    nz_fac: %" PetscInt_FMT "\n", info.nz_fac));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    fill: %.5e\n", info.fill));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    time: %.5e\n", info.time));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    factor_type: %s\n", info.factor_type));

  // setup local subdom mats using the still global is
  PetscCall(MatCreateSubMatrices(A, data->sz, data->is, data->is, MAT_INITIAL_MATRIX, &data->mat));

  // setup local rank data structures -> makes is local
  PetscCall(ISGetLocalSize(data->rank_is, &size));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, size, &data->rank_sol));
  PetscCall(VecScatterCreate(gtemp, data->rank_is, data->rank_sol, NULL, &data->rank_sc));

  // setup local subdom ksp and scatters
  if (data->print_local) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  local:\n"));
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(net2as_setup_ds(pc, PETSC_COMM_SELF, "", &data->ksp[i], &data->mat[i], &data->sol[i], &info));
    PetscCall(VecScatterCreate(data->rank_sol, data->local_is[i], data->sol[i], NULL, &data->sc[i]));
    t_loc += info.time;
    if (data->print_local) {
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "    - size: %" PetscInt_FMT "\n", info.size));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      nz_mat: %" PetscInt_FMT "\n", info.nz_mat));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      nz_fac: %" PetscInt_FMT "\n", info.nz_fac));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      fill: %.5e\n", info.fill));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      time: %.5e\n", info.time));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      factor_type: %s\n", info.factor_type));
    }
  }
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  PetscCallMPI(MPI_Reduce(&t_loc, &t_max, 1, MPIU_REAL, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&t_loc, &t_sum, 1, MPIU_REAL, MPI_SUM, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  load_bal:\n    est: %.5e\n    mes: %.5e\n", (double)data->bal_est, (double)t_max / t_sum * size));

  PetscCall(MatCOO_Free(&coo));
  PetscCall(MatDestroy(&coarse_basis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));

  PetscCall(VecScatterBegin(data->rank_sc, x, data->rank_sol, INSERT_VALUES, SCATTER_FORWARD));

  // coarse
  PetscCall(MatMultTranspose(data->cb, x, data->csol));
  PetscCall(KSPSolve(data->cksp, data->csol, data->csol));
  PetscCall(MatMult(data->cb, data->csol, y));

  PetscCall(VecScatterEnd(data->rank_sc, x, data->rank_sol, INSERT_VALUES, SCATTER_FORWARD));

  // local
  // NOTE: separate loops are required as otherwise the reverse scatter would modify the other rank entries
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(VecScatterBegin(data->sc[i], data->rank_sol, data->sol[i], INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(data->sc[i], data->rank_sol, data->sol[i], INSERT_VALUES, SCATTER_FORWARD));
  }
  PetscCall(VecZeroEntries(data->rank_sol));
  for (PetscInt i = 0; i < data->sz; i++)
    PetscCall(KSPSolve(data->ksp[i], data->sol[i], data->sol[i]));
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(VecScatterBegin(data->sc[i], data->sol[i], data->rank_sol, ADD_VALUES, SCATTER_REVERSE));
    PetscCall(VecScatterEnd(data->sc[i], data->sol[i], data->rank_sol, ADD_VALUES, SCATTER_REVERSE));
  }

  // rank -> global
  PetscCall(VecScatterBegin(data->rank_sc, data->rank_sol, y, ADD_VALUES, SCATTER_REVERSE));
  PetscCall(VecScatterEnd(data->rank_sc, data->rank_sol, y, ADD_VALUES, SCATTER_REVERSE));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
  if (!isascii) goto end;

  PetscCall(PetscViewerASCIIPrintf(viewer, "--- COARSE ---\n"));
  PetscCall(MatView(data->cb, viewer));
  PetscCall(MatView(data->cmat, viewer));
  PetscCall(KSPView(data->cksp, viewer));

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
  const char* part = "q1";
  const char* load = "gr";

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&data));
  pc->data = (void*)data;

  // minimal for Q1
  data->p[0] = data->p[1] = 1;
  data->eps = 1e-10;
  data->mult_bound = 10;
  data->delta = 2;
  data->pux_dim = 0;
  memcpy(data->cb_type, part, strlen(part)+1);
  memcpy(data->load_type, load, strlen(load)+1);

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

#if 0
PetscErrorCode PCNet2ASGetQuotientGraph(PC pc, Mat *q) {
  PC_Net2AS *data;
  Mat local_adj;
  PetscInt n_rows, n_cols;
  const PetscInt *ia, *ja;
  PetscBool done;
  PetscSF sf;
  PetscLayout col_layout;
  ISLocalToGlobalMapping col_lgmap;
  const PetscInt *col_globals, *part_idx;
  PetscInt *leaf_part;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  PetscAssertPointer(q, 2);
  data = (PC_Net2AS *)pc->data;

  // get local CSR
  PetscCall(MatGetLocalMat(data->adj, MAT_INITIAL_MATRIX, &local_adj));
  PetscCall(MatGetRowIJ(local_adj, 0, PETSC_FALSE, PETSC_TRUE, &n_rows, &ia, &ja, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "MatGetRowIJ failed");

  // get column local-to-global mapping and its size
  PetscCall(MatGetLocalToGlobalMapping(data->adj, NULL, &col_lgmap));
  PetscCall(ISLocalToGlobalMappingGetSize(col_lgmap, &n_cols));
  PetscCall(ISLocalToGlobalMappingGetIndices(col_lgmap, &col_globals));

  // build SF for ghost exchange
  PetscCall(MatGetLayouts(data->adj, &col_layout, NULL));
  PetscCall(PetscSFCreate(PETSC_COMM_WORLD, &sf));
  PetscCall(PetscSFSetGraphLayout(sf, col_layout, n_cols, NULL, PETSC_COPY_VALUES, col_globals));

  // scatter partition values to local columns
  PetscCall(ISGetIndices(data->partition, &part_idx));
  PetscCall(PetscMalloc1(n_cols, &leaf_part));
  PetscCall(PetscSFBcastBegin(sf, MPIU_INT, part_idx, leaf_part, MPI_REPLACE));
  PetscCall(PetscSFBcastEnd(sf, MPIU_INT, part_idx, leaf_part, MPI_REPLACE));

  // build quotient graph
  PetscCall(MatCreate(PETSC_COMM_WORLD, q));
  PetscCall(MatSetSizes(*q, PETSC_DECIDE, PETSC_DECIDE, data->n_parts, data->n_parts));
  PetscCall(MatSetUp(*q));

  for (PetscInt i = 0; i < n_rows; i++) {
    for (PetscInt k = ia[i]; k < ia[i + 1]; k++) {
      PetscInt j = ja[k];
      if (leaf_part[i] != leaf_part[j])
        PetscCall(MatSetValue(*q, leaf_part[i], leaf_part[j], 1.0, INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(*q, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*q, MAT_FINAL_ASSEMBLY));

  // cleanup
  PetscCall(ISRestoreIndices(data->partition, &part_idx));
  PetscCall(ISLocalToGlobalMappingRestoreIndices(col_lgmap, &col_globals));
  PetscCall(MatRestoreRowIJ(local_adj, 0, PETSC_FALSE, PETSC_TRUE, &n_rows, &ia, &ja, &done));
  PetscCall(MatDestroy(&local_adj));
  PetscCall(PetscFree(leaf_part));
  PetscCall(PetscSFDestroy(&sf));

  PetscFunctionReturn(PETSC_SUCCESS);
}


// 4-color a small planar graph (SeqAIJ on rank 0) using DSatur + backtracking.
// colors[] must be preallocated with size p, filled with result 0..3.
PetscErrorCode PCNet2ASColorPlanarGraph4(Mat Q, PetscInt p, PetscInt *colors) {
  const PetscInt *ia, *ja;
  PetscBool done;
  PetscInt *saturation, *avail; // avail: bitmask of available colors per vertex
  PetscInt colored = 0;

  PetscFunctionBegin;
  PetscCall(MatGetRowIJ(Q, 0, PETSC_FALSE, PETSC_TRUE, &p, &ia, &ja, &done));
  PetscCheck(done, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "MatGetRowIJ failed");

  PetscCall(PetscMalloc1(p, &saturation));
  PetscCall(PetscMalloc1(p, &avail));
  for (PetscInt i = 0; i < p; i++) {
    colors[i] = -1;
    saturation[i] = 0;
    avail[i] = 0xF; // bits 0..3 = colors 0..3 available
  }

  while (colored < p) {
    // DSatur: pick uncolored vertex with highest saturation, break ties by degree
    PetscInt best = -1, best_sat = -1, best_deg = -1;
    for (PetscInt i = 0; i < p; i++) {
      if (colors[i] >= 0) continue;
      PetscInt deg = ia[i + 1] - ia[i];
      if (saturation[i] > best_sat || (saturation[i] == best_sat && deg > best_deg)) {
        best = i;
        best_sat = saturation[i];
        best_deg = deg;
      }
    }

    // pick lowest available color (try 4, fall back to 5)
    PetscInt c = -1;
    for (PetscInt k = 0; k < 5; k++) {
      if (avail[best] & (1 << k)) { c = k; break; }
    }
    PetscCheck(c >= 0, PETSC_COMM_SELF, PETSC_ERR_PLIB,
      "5-coloring failed at vertex %" PetscInt_FMT, best);
    if (c == 4)
      PetscCall(PetscPrintf(PETSC_COMM_SELF,
        "WARNING: DSatur used 5th color at vertex %" PetscInt_FMT "\n", best));

    colors[best] = c;
    colored++;

    // update neighbors
    for (PetscInt k = ia[best]; k < ia[best + 1]; k++) {
      PetscInt nb = ja[k];
      if (colors[nb] < 0) {
        PetscInt old_avail = avail[nb];
        avail[nb] &= ~(1 << c);
        if (avail[nb] != old_avail) saturation[nb]++;
      }
    }
  }

  PetscCall(MatRestoreRowIJ(Q, 0, PETSC_FALSE, PETSC_TRUE, &p, &ia, &ja, &done));
  PetscCall(PetscFree(saturation));
  PetscCall(PetscFree(avail));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
