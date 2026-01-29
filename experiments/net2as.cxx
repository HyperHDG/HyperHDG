#include "net2as.hxx"
#include "prin2.hxx"
#include <petsc/private/pcimpl.h>
#include <petsc/private/hashmapi.h>
#include <petscviewerhdf5.h>

// PROBLEMS: if we have unequal number of local problems per rank, then we need to fill is sol and scatter with dummy/ empty stuff,
//   else we get a deadlock
// REFACTOR: use one rank scattering context, one rank is and one rank sol, then work with subvectors

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
  // number of local data structures (should be prod(p))
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
  // part_type one of "q1", "alg"
  char part_type[10];
  // load_type one of "rr", "gr"
  // rr - naive round robin load balancing
  // gr - simplest greedy load balancing
  char load_type[10];

  // network information

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
  Vec* sol;
  VecScatter* sc;
};

PetscErrorCode net2as_alloc_ds(PC_Net2AS *data, PetscInt sz) {
  PetscFunctionBegin;
  PetscCall(PetscMalloc5(sz, &data->ksp, sz, &data->is, sz,
    &data->sol, sz, &data->sc, sz, &data->sd_gids));
  data->sz = sz;
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

  // local subdomain resources
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(KSPDestroy(data->ksp+i));
    PetscCall(VecDestroy(data->sol+i));
    PetscCall(ISDestroy(data->is+i));
  }
  PetscCall(MatDestroySubMatrices(data->sz, &data->mat));
  PetscCall(PetscFree5(data->ksp, data->is, data->sol, data->sc, data->sd_gids));

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
  PetscCall(PetscOptionsString("-net2as_part_type", "subdomain partition type", NULL, data->part_type, data->part_type, sizeof(data->part_type), &set));
  PetscCall(PetscOptionsString("-net2as_load_type", "subdomain load balancing type", NULL, data->load_type, data->load_type, sizeof(data->load_type), &set));
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

  PetscCheck(is_size == n, comm, PETSC_ERR_ARG_SIZ, "types_points size %d != points size %d", is_size, n);
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
  PetscCall(MatSetSizes(data->adj, nn, nn, n, n));
  PetscCall(MatSetOptionsPrefix(data->adj, "adj_"));
  PetscCall(MatSetPreallocationCOO(data->adj, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(data->adj, coo.vals, INSERT_VALUES));
  PetscCall(MatCOO_Free(&coo));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_domain:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  points: %" PetscInt_FMT "\n", n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  edges: %" PetscInt_FMT "\n", m));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode net2as_setup_ds(PC pc, MPI_Comm comm, KSP *ksp, Mat *mat, Vec *sol, PetscLogDouble *t) {
  PC subpc;
  PetscLogDouble t0, t1;
  const char* prefix;

  PetscCall(PetscTime(&t0));
  PetscFunctionBegin;
  PetscCall(KSPCreate(comm, ksp));
  PetscCall(KSPSetType(*ksp, KSPPREONLY));
  PetscCall(KSPGetPC(*ksp, &subpc));
  PetscCall(PCSetType(subpc, PCCHOLESKY));

  PetscCall(PCGetOptionsPrefix(pc, &prefix));
  PetscCall(KSPSetOptionsPrefix(*ksp, prefix));
  PetscCall(KSPAppendOptionsPrefix(*ksp, "net2as_"));
  PetscCall(PCSetOptionsPrefix(subpc, prefix));
  PetscCall(PCAppendOptionsPrefix(subpc, "net2as_"));
  PetscCall(PCSetFromOptions(subpc));
  PetscCall(KSPSetFromOptions(*ksp));

  PetscCall(MatCreateVecs(*mat, sol, NULL));
  PetscCall(KSPSetOperators(*ksp, *mat, *mat));
  PetscCall(KSPSetUp(*ksp));
  PetscCall(PetscTime(&t1));

  *t = t1-t0;
  PetscFunctionReturn(0);
}

// simplest of all greedy load balancing strategies
PetscErrorCode net2as_loadbalance_greedy(MPI_Comm comm, PetscInt *weights, PetscInt *assignments, PetscInt count) {
  int size;
  PetscHeap loads;

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
    PetscCall(PetscHeapAdd(loads, r, load));
  }
  // PetscCall(PetscHeapView(loads, NULL));
  PetscCall(PetscHeapDestroy(&loads));
  PetscFunctionReturn(0);
}

// naive round robin scheduling
PetscErrorCode net2as_loadbalance_round_robin(MPI_Comm comm, PetscInt *weights, PetscInt *assignments, PetscInt count) {
  int size;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &size));
  for (PetscInt i = 0; i < count; i++) assignments[i] = i % size;
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
  PetscInt p = data->p[0]*data->p[1], sd_count, sd_total_size, off, start;
  PetscInt *sd2lcounts, *sd2gcounts, *sd2rank, *rank2scount, *rank2rcount, *coo2rank;
  MPI_Request *reqs;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCall(PetscCalloc7(p, &sd2lcounts, p, &sd2gcounts, p, &sd2rank,
    size, &rank2scount, size, &rank2rcount, cb->nnz, &coo2rank, 4*size, &reqs));

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
      PetscCall(net2as_loadbalance_round_robin(comm, sd2gcounts, sd2rank, p));
    else if (strcmp(data->load_type, "gr") == 0)
      PetscCall(net2as_loadbalance_greedy(comm, sd2gcounts, sd2rank, p));
    else
      PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
        "unexptected load balancing type '%s', expected one of 'rr', 'gr'", data->load_type);
  }
  // send this assignment to all ranks
  PetscCallMPI(MPI_Bcast(sd2rank, p, MPIU_INT, 0, comm));

  // PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "sd2rank", sd2rank, p));

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
  for (PetscInt i = 0; i < cb->nnz; i++) coo2rank[i] = sd2rank[cb->cols[i]];
  PetscCall(PetscSortIntWithArrayPair(cb->nnz, coo2rank, cb->rows, cb->cols));

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
  PetscCheck(off == sd_count, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "detected '%d' subdomains, expected '%d'", off, sd_count);

  PetscCall(PetscFree7(sd2lcounts, sd2gcounts, sd2rank,
    rank2rcount, rank2scount, coo2rank, reqs));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_is_blocked(PC_Net2AS *data) {
  PetscFunctionBegin;
  for (PetscInt i = 0; i < data->sz; i++) {
    const PetscInt *inds;
    PetscInt sz;
    IS is;
    PetscCall(ISGetIndices(data->is[i], &inds));
    PetscCall(ISGetLocalSize(data->is[i], &sz));
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, sz, inds, PETSC_COPY_VALUES, &is));
    PetscCall(ISRestoreIndices(data->is[i], &inds));
    PetscCall(ISDestroy(&data->is[i]));
    data->is[i] = is;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_rank_is(PC_Net2AS *data, MatCOO *sd) {
  PetscInt size = sd->nnz, *rank_is;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(sd->nnz, &rank_is));
  PetscCall(PetscArraycpy(rank_is, sd->cols, sd->nnz));
  PetscCall(PetscSortRemoveDupsInt(&size, rank_is));
  PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, size, rank_is, PETSC_COPY_VALUES, &data->rank_is));
  PetscCall(PetscFree(rank_is));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_q1(PC_Net2AS *data, MatCOO *coo, MatCOO *sd) {
  PetscReal h[2], min[2], max[2], eps = data->eps;
  PetscInt vstart, vend, size;
  const PetscInt *types;
  std::span<PetscReal> vspan;

  PetscFunctionBegin;
  for (PetscInt i = 0; i < 2; i++) {
    PetscCall(VecStrideMin(data->points, i, NULL, min+i));
    PetscCall(VecStrideMax(data->points, i, NULL, max+i));
    h[i] = (max[i]-min[i])/(data->p[i]+1);
  }

  PetscCall(ISGetIndices(data->types_points, &types));
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

    if (i>data->p[0] || j>data->p[1] || types[n] != 0) continue;
    if (i>0 && j>0 && (1-xx)*(1-yy) > eps)
      PetscCall(MatCOO_Push(coo, vstart+n, (j-1)*data->p[0]+(i-1), (1-xx)*(1-yy)));
    if (i < data->p[0] && j > 0 && PetscAbs(xx*(1-yy)) > eps)
      PetscCall(MatCOO_Push(coo, vstart+n, (j-1)*data->p[0]+i, xx*(1-yy)));
    if (i > 0 && j < data->p[1] && PetscAbs((1-xx)*yy) > eps)
      PetscCall(MatCOO_Push(coo, vstart+n, j*data->p[0]+(i-1), (1-xx)*yy));
    if (i < data->p[0] && j < data->p[1] && PetscAbs(xx*yy) > eps)
      PetscCall(MatCOO_Push(coo, vstart+n, j*data->p[0]+i, xx*yy));
  }
  PetscCall(VecRestoreSpan(data->points, vspan));

  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, coo, sd));
  PetscCall(net2as_make_is_blocked(data));
  PetscCall(net2as_make_rank_is(data, sd));

  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_alg(PC_Net2AS *data, MatCOO *coo, MatCOO *sd) {
  MatPartitioning p_ctx;
  IS partition;
  PetscInt p = data->p[0]*data->p[1], lsz_part, new_cap, vstart, vend;
  const PetscInt *inds, *types;
  PetscHMapI counts;

  PetscFunctionBegin;
  PetscCall(MatPartitioningCreate(PETSC_COMM_WORLD, &p_ctx));
  PetscCall(MatPartitioningSetAdjacency(p_ctx, data->adj));
  PetscCall(MatPartitioningSetNParts(p_ctx, p));
  PetscCall(MatPartitioningSetFromOptions(p_ctx));
  PetscCall(MatPartitioningApply(p_ctx, &partition));
  PetscCall(MatPartitioningDestroy(&p_ctx));

  PetscCall(ISGetIndices(partition, &inds));
  PetscCall(ISGetIndices(data->types_points, &types));
  PetscCall(ISGetLocalSize(partition, &lsz_part));
  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  vstart /= 3; vend /= 3;
  PetscCheck(vend-vstart == lsz_part, PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
    "parallel layout of data->points must match that of the partition, data->points size '%d', "
    "partition size '%d'", vend-vstart, lsz_part);
  PetscCall(MatCOO_Alloc(coo, lsz_part));
  for (PetscInt i = 0; i < lsz_part; i++)
    if (types[i] == 0) PetscCall(MatCOO_Push(coo, vstart+i, inds[i], 1.));
  PetscCall(ISRestoreIndices(partition, &inds));
  PetscCall(ISRestoreIndices(data->types_points, &types));

  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, coo, sd));
  PetscCall(MatIncreaseOverlap(data->adj, data->sz, data->is, data->delta));
  PetscCall(net2as_make_is_blocked(data));
  PetscCall(net2as_make_rank_is(data, sd));

  new_cap = 0;
  PetscCall(PetscHMapICreate(&counts));
  for (PetscInt s = 0; s < data->sz; s++) {
    PetscInt sz, val;
    PetscCall(ISGetIndices(data->is[s], &inds)); // these will be global indices, not local ones which we need
    PetscCall(ISGetLocalSize(data->is[s], &sz));
    for (PetscInt i = 0; i < sz; i++) {
      PetscCall(PetscHMapIGetWithDefault(counts, inds[i], 0, &val));
      PetscCall(PetscHMapISet(counts, inds[i], val + 1));
    }
    PetscCall(ISRestoreIndices(data->is[s], &inds));
    new_cap += sz;
  }

  PetscCall(MatCOO_Free(coo));
  PetscCall(MatCOO_Alloc(coo, new_cap));

  for (PetscInt s = 0; s < data->sz; s++) {
    PetscInt sz, count;
    PetscCall(ISGetIndices(data->is[s], &inds));
    PetscCall(ISGetLocalSize(data->is[s], &sz));
    for (PetscInt i = 0; i < sz; i++) {
      PetscCall(PetscHMapIGet(counts, inds[i], &count));
      PetscCall(MatCOO_Push(coo, inds[i], data->sd_gids[i], 1./count));
    }
    PetscCall(ISRestoreIndices(data->is[s], &inds));
  }

  PetscCall(PetscHMapIDestroy(&counts));

  PetscFunctionReturn(0);
}

PetscErrorCode net2as_make_is_local(PC_Net2AS *data, MatCOO *sd) {
  ISLocalToGlobalMapping l2g;

  PetscFunctionBegin;
  PetscCall(ISLocalToGlobalMappingCreateIS(data->rank_is, &l2g));
  for (PetscInt i = 0; i < data->sz; i++) {
    const PetscInt *global;
    PetscInt *local;
    PetscInt n_global, n_local;
    IS is;

    PetscCall(ISBlockGetLocalSize(data->is[i], &n_global));
    PetscCall(ISBlockGetIndices(data->is[i], &global));

    PetscCall(PetscMalloc1(n_global, &local));
    PetscCall(ISGlobalToLocalMappingApply(l2g, IS_GTOLM_DROP, n_global, global, &n_local, local));
    PetscCheck(n_global == n_local, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
      "rank_is local to global mapping inds dropped: expected %" PetscInt_FMT ", got %" PetscInt_FMT, n_global, n_local);
    PetscCall(ISBlockRestoreIndices(data->is[i], &global));
    PetscCall(ISCreateBlock(PETSC_COMM_SELF, data->bs, n_global, local, PETSC_OWN_POINTER, &is));
    PetscCall(ISDestroy(&data->is[i]));
    data->is[i] = is;
  }
  PetscCall(ISLocalToGlobalMappingDestroy(&l2g));

  PetscFunctionReturn(0);
}

PetscErrorCode PCSetup_Net2AS(PC pc) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  Vec gtemp;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscInt size, msize, n_cols = data->p[0] * data->p[1], n, m;
  Mat coarse_basis, A;
  MatType type;
  int comm_size;
  MatCOO coo, sd;
  PetscLogDouble t;

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

  if (strcmp(data->part_type, "q1") == 0)
    PetscCall(net2as_cb_q1(data, &coo, &sd));
  else if (strcmp(data->part_type, "alg") == 0)
    PetscCall(net2as_cb_alg(data, &coo, &sd));
  else
    PetscCheck(false, PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "unsupported type '%s', muse be one of 'q1', 'alg'", data->part_type);

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
  PetscCall(MatCreateMAIJ(coarse_basis, data->bs, &data->cb)); // expanded by block size

  // setup coarse global data structures
  PetscCall(MatPtAP(A, data->cb, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->cmat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  part_type: %s\n", data->part_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  load_type: %s\n", data->load_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  bs: %" PetscInt_FMT "\n", data->bs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  p: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", data->p[0], data->p[1]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  sz: %" PetscInt_FMT "\n", n_cols));
  PetscCall(MatGetSize(A, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  global_size: %" PetscInt_FMT "\n", m));
  PetscCall(MatGetSize(data->cmat, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  coarse:\n    size: %" PetscInt_FMT "\n", m));
  PetscCall(net2as_setup_ds(pc, PETSC_COMM_WORLD, &data->cksp, &data->cmat, &data->csol, &t));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    time: %.5e\n", t));

  // setup local subdom mats using the still global is
  PetscCall(MatCreateSubMatrices(A, data->sz, data->is, data->is, MAT_INITIAL_MATRIX, &data->mat));

  // setup local rank data structures -> makes is local
  PetscCall(net2as_make_is_local(data, &sd));
  PetscCall(ISGetLocalSize(data->rank_is, &size));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, size, &data->rank_sol));
  PetscCall(VecScatterCreate(gtemp, data->rank_is, data->rank_sol, NULL, &data->rank_sc));

  // setup local subdom ksp and scatters
  if (data->print_local) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  local:\n"));
  for (PetscInt i = 0; i < data->sz; i++) {
    PetscCall(net2as_setup_ds(pc, PETSC_COMM_SELF, &data->ksp[i], &data->mat[i], &data->sol[i], &t));
    PetscCall(VecScatterCreate(data->rank_sol, data->is[i], data->sol[i], NULL, &data->sc[i]));
    if (data->print_local) {
      PetscInt local_size;
      PetscCall(ISGetLocalSize(data->is[i], &local_size));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "    - size: %" PetscInt_FMT "\n", local_size));
      PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      time: %.5e\n", t));
    }
  }
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

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
  memcpy(data->part_type, part, strlen(part)+1);
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
