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
  // overlap target as a fraction of the (weighted) graph diameter; delta = overlap_frac * diam
  PetscReal overlap_frac;
  // absolute overlap distance override in edge-length units; if > 0, used instead of overlap_frac*diam
  PetscReal overlap_abs;
  // the number of dimension to extend the PU to (linear enrichment x, y[, z])
  PetscInt pux_dim;
  // (cb_pu only) add the bilinear cross term xy to the coarse space, i.e. reproduce {1,x,y,xy} like
  // q1 rather than only the linear {1,x,y}; not the full quadratic (no x^2, y^2)
  PetscBool pu_xy;
  // number of coarse basis functions per subdomain (1 constant + pux_dim linear + pu_xy cross)
  PetscInt cb_ncomp;
  // cb_type one of "q1", "pu"
  char cb_type[10];
  // load_type one of "rr", "gr"
  // rr - naive round robin load balancing
  // gr - simplest greedy load balancing
  char load_type[10];
  // apply no coarse correction
  PetscBool nocoarse;
  // (cb_q1 only) trim the coarse DoFs that peak on the domain boundary
  PetscBool cb_trim;
  // (non-overlapping) subdomain diameter
  PetscReal H;

  // network information

  // path to domain file
  char domain[PATH_MAX];
  // when true, points+adj were injected via PCNet2ASSetDomain (already redistributed to match the
  // system matrix) and the domain file is not read.
  PetscBool domain_set;
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
  PetscCall(PetscOptionsReal("-net2as_overlap_frac", "overlap distance as a fraction of the weighted graph diameter", NULL, data->overlap_frac, &data->overlap_frac, &set));
  PetscCall(PetscOptionsReal("-net2as_overlap_abs", "absolute overlap distance in edge-length units; overrides overlap_frac when > 0", NULL, data->overlap_abs, &data->overlap_abs, &set));
  PetscCall(PetscOptionsInt("-net2as_pux_dim", "number of dimension to extend the pu by", NULL, data->pux_dim, &data->pux_dim, &set));
  PetscCall(PetscOptionsBool("-net2as_pu_xy", "add the bilinear cross term xy to the pu coarse space ({1,x,y,xy} like q1)", NULL, data->pu_xy, &data->pu_xy, &set));
  PetscCall(PetscOptionsString("-net2as_cb_type", "subdomain partition type", NULL, data->cb_type, data->cb_type, sizeof(data->cb_type), &set));
  PetscCall(PetscOptionsString("-net2as_load_type", "subdomain load balancing type", NULL, data->load_type, data->load_type, sizeof(data->load_type), &set));
  PetscCall(PetscOptionsBool("-net2as_nocoarse", "apply no coarse correction", NULL, data->nocoarse, &data->nocoarse, &set));
  PetscCall(PetscOptionsBool("-net2as_cb_trim", "trim the cb_q1 coarse DoFs that peak on the domain boundary (coarse space only, BC-conforming; the subdomain cover keeps the boundary patches)", NULL, data->cb_trim, &data->cb_trim, &set));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscBool net2as_is_dirichlet(PetscInt type, PetscBool wave) {
  if (type == 0) return PETSC_FALSE;
  if (wave && (type & (1u << 6))) return PETSC_FALSE;  // static-only, free in wave
  return PETSC_TRUE;
}

// Inject an already-redistributed domain (points + edges in the partition's global numbering) so
// that net2as's points/adjacency conform to the system matrix's parallel layout. Replaces reading
// the (un-partitioned) domain file. coords: n_owned_nodes*sdim, row-major owned points in global
// order; edges_global: 2*n_owned_edges global hypernode index pairs (this rank's owned edges).
PetscErrorCode PCNet2ASSetDomain(PC pc,
                                 PetscInt n_owned_nodes, PetscInt n_global_nodes, PetscInt sdim,
                                 const PetscReal *coords,
                                 PetscInt n_owned_edges, const PetscInt *edges_global) {
  PC_Net2AS *data = (PC_Net2AS*)pc->data;
  MPI_Comm comm = PetscObjectComm((PetscObject)pc);
  PetscScalar *arr;
  MatCOO coo;

  PetscFunctionBeginUser;
  PetscCall(VecCreate(comm, &data->points));
  PetscCall(PetscObjectSetName((PetscObject)data->points, "points"));
  PetscCall(VecSetSizes(data->points, n_owned_nodes * sdim, n_global_nodes * sdim));
  PetscCall(VecSetBlockSize(data->points, sdim));
  PetscCall(VecSetType(data->points, VECMPI));
  PetscCall(VecGetArray(data->points, &arr));
  for (PetscInt i = 0; i < n_owned_nodes * sdim; i++) arr[i] = coords[i];
  PetscCall(VecRestoreArray(data->points, &arr));

  PetscCall(MatCOO_Alloc(&coo, 2 * n_owned_edges));
  for (PetscInt i = 0; i < n_owned_edges; i++) {
    PetscCall(MatCOO_Push(&coo, edges_global[2*i], edges_global[2*i+1], 1.));
    PetscCall(MatCOO_Push(&coo, edges_global[2*i+1], edges_global[2*i], 1.));
  }
  PetscCall(MatCreate(comm, &data->adj));
  PetscCall(MatSetType(data->adj, MATMPIAIJ));
  PetscCall(MatSetSizes(data->adj, n_owned_nodes, n_owned_nodes, n_global_nodes, n_global_nodes));
  PetscCall(MatSetOptionsPrefix(data->adj, "net2as_adj_"));
  PetscCall(MatSetPreallocationCOO(data->adj, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(data->adj, coo.vals, INSERT_VALUES));
  PetscCall(MatCOO_Free(&coo));

  data->domain_set = PETSC_TRUE;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_domain:\n  injected: true\n  points: %"
                        PetscInt_FMT "\n", n_global_nodes));
  PetscFunctionReturn(PETSC_SUCCESS);
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
    PetscCallMPI(MPI_Isend(tmp_rows+off, rank2scount[r], MPIU_INT, r, tag_vid, comm, &reqs[2*(size+r)]));
    PetscCallMPI(MPI_Isend(tmp_cols+off, rank2scount[r], MPIU_INT, r, tag_sid, comm, &reqs[2*(size+r)+1]));
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
  data->H = PetscMax(h[0], h[1]);
  data->overlap_abs = data->H;
  data->overlap_frac = 1;

  // Optionally trim the coarse DoFs that peak on the domain boundary: the outer ring of the
  // tensor-product Q1 grid (i in {0, ns[0]-1} or j in {0, ns[1]-1}). These basis functions peak on
  // the Dirichlet boundary where the solution is fixed; dropping them makes the coarse space
  // conform to the homogeneous boundary condition (gortz.pdf's V_0: for an aligned coarse grid the
  // remaining interior hats vanish on the boundary). The trim applies to the COARSE SPACE only —
  // the subdomain cover below keeps the boundary patches, matching the paper's decomposition
  // (all partition-of-unity patches, coarse space in V).
  // col_remap maps each old coarse column to its compacted index, or -1 if trimmed.
  PetscInt *col_remap = NULL;
  PetscInt n_coarse_kept = n_coarse;
  if (data->cb_trim) {
    n_coarse_kept = 0;
    PetscCall(PetscMalloc1(n_coarse, &col_remap));
    for (PetscInt j = 0; j < ns[1]; j++)
      for (PetscInt i = 0; i < ns[0]; i++) {
        PetscBool bdry = (PetscBool)(i == 0 || i == ns[0]-1 || j == 0 || j == ns[1]-1);
        col_remap[j*ns[0]+i] = bdry ? -1 : n_coarse_kept++;
      }
  }

  PetscCall(VecGetOwnershipRange(data->points, &vstart, &vend));
  n_local = (vend - vstart) / 3;
  vstart /= 3;

  PetscCall(VecGetArray(data->points, &points));
  data->n_coarse = n_coarse;
  data->cb_ncomp = 1; // q1 has no per-subdomain enrichment

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
      // Skip exact-zero weights: they would put nodes lying exactly on a patch boundary into
      // the far-side subdomains (closed instead of open hat supports). On coarse meshes aligned
      // with a regular grid network that couples same-color patches and raises lambda_max of the
      // preconditioned operator from 4 to 6-9; gortz.pdf's Table 2 grid rates are only
      // reproduced with open supports.
      if (pts[l].w == 0) continue;
      PetscCall(MatCOO_Push(coo, row, col, pts[l].w));
      // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%.2e, %.2e, %d, %d, | %d, %d, %.2e\n", x, y, i, j , row, col, pts[l].w));
    }
  }

  PetscCall(VecRestoreArray(data->points, &points));

  // subdomains = supports of ALL hats (including the boundary-peaked ones), so the cover has no
  // thin spots along the boundary even when the coarse space is trimmed
  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, coo, sd));
  PetscCall(net2as_make_rank_is(data->is, data->sz, data->bs, &data->rank_is));
  PetscCall(net2as_make_is_local(data, sd));
  PetscCall(net2as_make_is_blocked(data));

  // cb_trim: now restrict the coarse basis to the kept interior hats
  if (data->cb_trim) {
    PetscInt nnz = 0;
    for (PetscInt t = 0; t < coo->nnz; t++) {
      PetscInt col = col_remap[coo->cols[t]];
      if (col < 0) continue; // trimmed boundary coarse DoF
      coo->rows[nnz] = coo->rows[t];
      coo->cols[nnz] = col;
      coo->vals[nnz] = coo->vals[t];
      nnz++;
    }
    coo->nnz = nnz;
    data->n_coarse = n_coarse_kept;
  }

  PetscCall(PetscFree(col_remap));
  PetscFunctionReturn(0);
}

// A distance large enough to mark "unreached" without overflowing under +edge_length additions.
#define NET2AS_BIG (PETSC_MAX_REAL / 4)

// Weighted graph view of data->adj for shortest-path computations. Edge weights are the physical
// edge lengths ||p_i - p_j|| derived from the node coordinates. Built on the rank-local CSR
// (MatGetLocalMat) whose columns are owned nodes [0,n_rows) followed by ghosts; a PetscSF over the
// adjacency column layout exchanges per-node values (coords, distance labels) with the owners.
struct Net2AS_WGraph {
  PetscInt   *ia, *ja;      // combined local CSR; column indices in [0,n_cols): owned then ghost
  PetscInt    n_rows;       // owned nodes
  PetscInt    n_cols;       // owned (n_rows) + ghost columns referenced locally
  PetscInt    vstart;       // global id of first owned node
  PetscInt    nglobal;      // global node count (iteration safeguard)
  PetscReal  *w;            // edge weights (lengths), length ia[n_rows]
  PetscReal  *coords;       // n_cols*3 coords of every local column (owned + ghost)
  PetscInt   *col_globals;  // n_cols, global node id of each local column
  PetscSF     sf;           // roots: adj column layout (node layout); leaves: local columns
};

PetscErrorCode net2as_wgraph_create(Mat adj, Vec points, Net2AS_WGraph *g) {
  Mat Ad, Ao;
  const PetscInt *garray, *iad, *jad, *iao, *jao;
  PetscInt n_local, n_ghost, nd, no, nnz, off = 0;
  PetscBool done;
  PetscLayout col_layout;
  const PetscScalar *parr;
  PetscReal *rootc;
  MPI_Datatype unit3;

  PetscFunctionBegin;
  PetscCall(MatGetOwnershipRange(adj, &g->vstart, NULL));
  PetscCall(MatGetSize(adj, &g->nglobal, NULL));

  // diagonal (owned columns) + off-diagonal (ghost columns) blocks, no copy
  PetscCall(MatMPIAIJGetSeqAIJ(adj, &Ad, &Ao, &garray));
  PetscCall(MatGetRowIJ(Ad, 0, PETSC_FALSE, PETSC_FALSE, &nd, &iad, &jad, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "MatGetRowIJ(diag) failed");
  PetscCall(MatGetRowIJ(Ao, 0, PETSC_FALSE, PETSC_FALSE, &no, &iao, &jao, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "MatGetRowIJ(offdiag) failed");
  n_local = nd;
  PetscCall(MatGetSize(Ao, NULL, &n_ghost));
  g->n_rows = n_local;
  g->n_cols = n_local + n_ghost;

  // local column -> global node id: owned columns [0,n_local), then ghosts via garray
  PetscCall(PetscMalloc1(g->n_cols, &g->col_globals));
  for (PetscInt c = 0; c < n_local; c++) g->col_globals[c] = g->vstart + c;
  for (PetscInt c = 0; c < n_ghost; c++) g->col_globals[n_local + c] = garray[c];

  // merge the two blocks into one CSR in the unified local column numbering
  nnz = iad[n_local] + iao[n_local];
  PetscCall(PetscMalloc2(n_local + 1, &g->ia, nnz, &g->ja));
  g->ia[0] = 0;
  for (PetscInt i = 0; i < n_local; i++) {
    for (PetscInt k = iad[i]; k < iad[i+1]; k++) g->ja[off++] = jad[k];
    for (PetscInt k = iao[i]; k < iao[i+1]; k++) g->ja[off++] = n_local + jao[k];
    g->ia[i+1] = off;
  }
  PetscCall(MatRestoreRowIJ(Ad, 0, PETSC_FALSE, PETSC_FALSE, &nd, &iad, &jad, &done));
  PetscCall(MatRestoreRowIJ(Ao, 0, PETSC_FALSE, PETSC_FALSE, &no, &iao, &jao, &done));

  // SF mapping each local column (leaf) to the owner of its global node (root, in adj's col layout)
  PetscCall(MatGetLayouts(adj, NULL, &col_layout));
  PetscCall(PetscSFCreate(PETSC_COMM_WORLD, &g->sf));
  PetscCall(PetscSFSetGraphLayout(g->sf, col_layout, g->n_cols, NULL, PETSC_COPY_VALUES, g->col_globals));

  // gather coordinates of owned + ghost columns
  MPI_Type_contiguous(3, MPIU_REAL, &unit3);
  MPI_Type_commit(&unit3);
  PetscCall(PetscMalloc1(n_local * 3, &rootc));
  PetscCall(VecGetArrayRead(points, &parr));
  for (PetscInt i = 0; i < n_local * 3; i++) rootc[i] = PetscRealPart(parr[i]);
  PetscCall(VecRestoreArrayRead(points, &parr));
  PetscCall(PetscMalloc1(g->n_cols * 3, &g->coords));
  PetscCall(PetscSFBcastBegin(g->sf, unit3, rootc, g->coords, MPI_REPLACE));
  PetscCall(PetscSFBcastEnd(g->sf, unit3, rootc, g->coords, MPI_REPLACE));
  PetscCall(PetscFree(rootc));
  MPI_Type_free(&unit3);

  // edge weights = euclidean lengths
  PetscCall(PetscMalloc1(g->ia[g->n_rows], &g->w));
  for (PetscInt i = 0; i < g->n_rows; i++) {
    for (PetscInt k = g->ia[i]; k < g->ia[i+1]; k++) {
      PetscInt j = g->ja[k];
      PetscReal d2 = 0;
      for (PetscInt c = 0; c < 3; c++) {
        PetscReal dc = g->coords[3*i+c] - g->coords[3*j+c];
        d2 += dc * dc;
      }
      g->w[k] = PetscSqrtReal(d2);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_wgraph_destroy(Net2AS_WGraph *g) {
  PetscFunctionBegin;
  PetscCall(PetscSFDestroy(&g->sf));
  PetscCall(PetscFree2(g->ia, g->ja));
  PetscCall(PetscFree(g->w));
  PetscCall(PetscFree(g->coords));
  PetscCall(PetscFree(g->col_globals));
  PetscFunctionReturn(0);
}

// Bounded multi-source weighted shortest path (Jacobi Bellman-Ford with ghost exchange).
// nsrc parallel distance fields are stored interleaved: Dloc[i*nsrc + s] is owned node i's distance
// to source set s. Caller initializes sources to 0 and all else to NET2AS_BIG. Candidate distances
// above delta_cap are not relaxed (pass PETSC_MAX_REAL for an uncapped solve).
PetscErrorCode net2as_sssp(Net2AS_WGraph *g, PetscInt nsrc, PetscReal delta_cap, PetscReal *Dloc) {
  PetscReal *leafD;
  MPI_Datatype unitn;
  PetscInt iter = 0;

  PetscFunctionBegin;
  MPI_Type_contiguous(nsrc, MPIU_REAL, &unitn);
  MPI_Type_commit(&unitn);
  PetscCall(PetscMalloc1(g->n_cols * nsrc, &leafD));

  while (1) {
    PetscInt changed = 0, gchanged;

    PetscCall(PetscSFBcastBegin(g->sf, unitn, Dloc, leafD, MPI_REPLACE));
    PetscCall(PetscSFBcastEnd(g->sf, unitn, Dloc, leafD, MPI_REPLACE));
    for (PetscInt i = 0; i < g->n_rows; i++) {
      for (PetscInt k = g->ia[i]; k < g->ia[i+1]; k++) {
        PetscInt j = g->ja[k];
        PetscReal wij = g->w[k];
        for (PetscInt s = 0; s < nsrc; s++) {
          PetscReal cand = leafD[j*nsrc + s] + wij;
          if (cand < Dloc[i*nsrc + s] && cand <= delta_cap) {
            Dloc[i*nsrc + s] = cand;
            changed = 1;
          }
        }
      }
    }
    PetscCallMPI(MPI_Allreduce(&changed, &gchanged, 1, MPIU_INT, MPI_LOR, PETSC_COMM_WORLD));
    if (!gchanged) break;
    if (++iter > g->nglobal + 2) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "WARNING: net2as_sssp hit iteration safeguard\n"));
      break;
    }
  }

  PetscCall(PetscFree(leafD));
  MPI_Type_free(&unitn);
  PetscFunctionReturn(0);
}

// Single-field weighted shortest path restricted to within-core edges: distances propagate only
// between nodes of the same partition (core). rowpart[i] = partition of owned node i; leafpart[c] =
// partition of local column c (owned + ghost, filled via SF). Caller seeds D (0 at sources, BIG else).
PetscErrorCode net2as_sssp_core(Net2AS_WGraph *g, const PetscInt *rowpart, const PetscInt *leafpart, PetscReal *D) {
  PetscReal *leafD;
  PetscInt iter = 0;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(g->n_cols, &leafD));
  while (1) {
    PetscInt changed = 0, gchanged;

    PetscCall(PetscSFBcastBegin(g->sf, MPIU_REAL, D, leafD, MPI_REPLACE));
    PetscCall(PetscSFBcastEnd(g->sf, MPIU_REAL, D, leafD, MPI_REPLACE));
    for (PetscInt i = 0; i < g->n_rows; i++) {
      PetscInt si = rowpart[i];
      for (PetscInt k = g->ia[i]; k < g->ia[i+1]; k++) {
        PetscInt j = g->ja[k];
        if (leafpart[j] != si) continue;          // stay inside the core
        PetscReal cand = leafD[j] + g->w[k];
        if (cand < D[i]) { D[i] = cand; changed = 1; }
      }
    }
    PetscCallMPI(MPI_Allreduce(&changed, &gchanged, 1, MPIU_INT, MPI_LOR, PETSC_COMM_WORLD));
    if (!gchanged) break;
    if (++iter > g->nglobal + 2) { PetscCall(PetscPrintf(PETSC_COMM_WORLD, "WARNING: net2as_sssp_core safeguard\n")); break; }
  }
  PetscCall(PetscFree(leafD));
  PetscFunctionReturn(0);
}

// Per-core weighted diameter via a double-sweep restricted to within-core edges: seed each core at
// its lowest-global-id node -> farthest node b_s -> diameter_s = max distance from b_s. diam[p] filled
// (>= 0; 0 for a singleton/empty core). rowpart = owned partition, leafpart = column partition (SF).
PetscErrorCode net2as_subdomain_diameters(Net2AS_WGraph *g, const PetscInt *rowpart,
                                          const PetscInt *leafpart, PetscInt p, PetscReal *diam) {
  PetscReal *D, *lmaxd, *gmaxd;
  PetscInt *seed, *lseed, *bnode, *lbnode;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(g->n_rows, &D));
  PetscCall(PetscMalloc4(p, &seed, p, &lseed, p, &bnode, p, &lbnode));
  PetscCall(PetscMalloc2(p, &lmaxd, p, &gmaxd));

  // seed_s = min global id in core s
  for (PetscInt s = 0; s < p; s++) lseed[s] = PETSC_MAX_INT;
  for (PetscInt i = 0; i < g->n_rows; i++) {
    PetscInt s = rowpart[i], gid = g->vstart + i;
    if (gid < lseed[s]) lseed[s] = gid;
  }
  PetscCallMPI(MPI_Allreduce(lseed, seed, p, MPIU_INT, MPI_MIN, PETSC_COMM_WORLD));

  // sweep 1 from seeds, then per-core argmax (max dist, min gid tie-break) = farthest node b_s
  for (PetscInt i = 0; i < g->n_rows; i++) D[i] = (g->vstart + i == seed[rowpart[i]]) ? 0 : NET2AS_BIG;
  PetscCall(net2as_sssp_core(g, rowpart, leafpart, D));
  for (PetscInt s = 0; s < p; s++) lmaxd[s] = -1;
  for (PetscInt i = 0; i < g->n_rows; i++) {
    PetscInt s = rowpart[i];
    if (D[i] < NET2AS_BIG/2 && D[i] > lmaxd[s]) lmaxd[s] = D[i];
  }
  PetscCallMPI(MPI_Allreduce(lmaxd, gmaxd, p, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
  for (PetscInt s = 0; s < p; s++) lbnode[s] = PETSC_MAX_INT;
  for (PetscInt i = 0; i < g->n_rows; i++) {
    PetscInt s = rowpart[i];
    if (D[i] < NET2AS_BIG/2 && D[i] == gmaxd[s] && g->vstart + i < lbnode[s]) lbnode[s] = g->vstart + i;
  }
  PetscCallMPI(MPI_Allreduce(lbnode, bnode, p, MPIU_INT, MPI_MIN, PETSC_COMM_WORLD));

  // sweep 2 from b_s, diameter_s = max distance
  for (PetscInt i = 0; i < g->n_rows; i++) D[i] = (g->vstart + i == bnode[rowpart[i]]) ? 0 : NET2AS_BIG;
  PetscCall(net2as_sssp_core(g, rowpart, leafpart, D));
  for (PetscInt s = 0; s < p; s++) lmaxd[s] = 0;
  for (PetscInt i = 0; i < g->n_rows; i++) {
    PetscInt s = rowpart[i];
    if (D[i] < NET2AS_BIG/2 && D[i] > lmaxd[s]) lmaxd[s] = D[i];
  }
  PetscCallMPI(MPI_Allreduce(lmaxd, diam, p, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));

  PetscCall(PetscFree(D));
  PetscCall(PetscFree4(seed, lseed, bnode, lbnode));
  PetscCall(PetscFree2(lmaxd, gmaxd));
  PetscFunctionReturn(0);
}

PetscErrorCode net2as_cb_pu(PC_Net2AS *data, MatCOO *coo, MatCOO *sd) {
  MatPartitioning p_ctx;
  MatPartitioningType p_type;
  char ptype[64];
  IS partition;
  PetscInt p = data->p[0]*data->p[1], lsz_part, vstart, vend, cut, nmemb, ncb;
  const PetscInt *part_owned;
  PetscInt *leafpart;
  Net2AS_WGraph g;
  PetscReal *diam_s, *delta_s, cap = 0, *Dloc, *cen = NULL, *ext = NULL;
  PetscReal dmin = PETSC_MAX_REAL, dmax = 0, dsum = 0;
  MatCOO memb;

  PetscFunctionBegin;

  data->n_coarse = p;
  data->cb_ncomp = data->pux_dim + 1 + (data->pu_xy ? 1 : 0);
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
  PetscCall(PetscStrncpy(ptype, p_type, sizeof(ptype))); // p_type points into p_ctx; copy before destroy
  PetscCall(MatPartitioningDestroy(&p_ctx));
  data->partition = partition; // kept for PCDestroy_Net2AS

  // build the weighted graph (edge lengths from node coordinates)
  PetscCall(net2as_wgraph_create(data->adj, data->points, &g));

  PetscCall(ISGetLocalSize(partition, &lsz_part));
  PetscCheck(vend-vstart == lsz_part, PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
    "parallel layout of data->points must match that of the partition, data->points size '%d', "
    "partition size '%d'", vend-vstart, lsz_part);
  PetscCheck(g.n_rows == lsz_part, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
    "weighted graph owns %" PetscInt_FMT " rows, partition owns %" PetscInt_FMT, g.n_rows, lsz_part);
  PetscCall(ISGetIndices(partition, &part_owned));

  // partition of each local column (owned + ghost), used for within-core distance computations
  PetscCall(PetscMalloc1(g.n_cols, &leafpart));
  PetscCall(PetscSFBcastBegin(g.sf, MPIU_INT, part_owned, leafpart, MPI_REPLACE));
  PetscCall(PetscSFBcastEnd(g.sf, MPIU_INT, part_owned, leafpart, MPI_REPLACE));

  // overlap distance is a fixed fraction of EACH subdomain's own weighted graph diameter (so the
  // condition-number bound ~1+H/delta stays constant as the subdomain count changes), not of the
  // global diameter. cap = max delta_s bounds the (shared) overlap SSSP.
  PetscCall(PetscMalloc2(p, &diam_s, p, &delta_s));
  PetscCall(net2as_subdomain_diameters(&g, part_owned, leafpart, p, diam_s));
  for (PetscInt s = 0; s < p; s++) {
    delta_s[s] = data->overlap_abs > 0 ? data->overlap_abs : data->overlap_frac * diam_s[s];
    if (delta_s[s] <= 0) delta_s[s] = PETSC_SMALL; // degenerate (singleton) core: core only, no overlap
    cap = PetscMax(cap, delta_s[s]);
    dmin = PetscMin(dmin, diam_s[s]); dmax = PetscMax(dmax, diam_s[s]); dsum += diam_s[s];
  }

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as_cb_pu:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  cut: %" PetscInt_FMT "\n", cut));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  part_type: %s\n", ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  pux_dim: %" PetscInt_FMT "\n", data->pux_dim));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  pu_xy: %s\n", data->pu_xy ? "true" : "false"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  cb_ncomp: %" PetscInt_FMT "\n", data->cb_ncomp));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  overlap_frac: %.5e\n", (double)data->overlap_frac));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  diam_subdomain: {min: %.5e, mean: %.5e, max: %.5e}\n",
                        (double)dmin, (double)(dsum / p), (double)dmax));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  delta_max: %.5e\n", (double)cap));

  data->H = dsum/p;

  // weighted shortest path from every node to each subdomain core, capped at the largest delta_s.
  // Dloc[i*p+s] = owned node i's distance to core s; node i is in subdomain s iff <= delta_s[s].
  PetscCall(PetscMalloc1(lsz_part * p, &Dloc));
  for (PetscInt i = 0; i < lsz_part * p; i++) Dloc[i] = NET2AS_BIG;
  for (PetscInt i = 0; i < lsz_part; i++) Dloc[i*p + part_owned[i]] = 0;
  PetscCall(net2as_sssp(&g, p, cap, Dloc));
  PetscCall(ISRestoreIndices(partition, &part_owned));

  // membership matrix (owned vertex, subdomain) defines the overlapping subdomains; distribute it
  nmemb = 0;
  for (PetscInt i = 0; i < lsz_part; i++)
    for (PetscInt s = 0; s < p; s++)
      if (Dloc[i*p + s] <= delta_s[s]) nmemb++;
  PetscCall(MatCOO_Alloc(&memb, nmemb));
  for (PetscInt i = 0; i < lsz_part; i++)
    for (PetscInt s = 0; s < p; s++)
      if (Dloc[i*p + s] <= delta_s[s]) PetscCall(MatCOO_Push(&memb, vstart+i, s, 1.));

  PetscCall(net2as_distribute_subdomains(PETSC_COMM_WORLD, data, &memb, sd));
  PetscCall(net2as_make_rank_is(data->is, data->sz, data->bs, &data->rank_is));
  PetscCall(net2as_make_is_local(data, sd));
  PetscCall(net2as_make_is_blocked(data));
  PetscCall(MatCOO_Free(&memb));

  // optional linear enrichment needs each subdomain's centroid and axis extent over its overlapping
  // support; reduce them globally since the support spans ranks.
  if (data->pux_dim > 0 || data->pu_xy) {
    PetscReal *csum, *cmin, *cmax, *gsum, *gmin, *gmax;
    PetscInt  *ccnt, *gcnt;
    PetscCall(PetscCalloc4(p*3, &csum, p, &ccnt, p*3, &cmin, p*3, &cmax));
    PetscCall(PetscMalloc4(p*3, &gsum, p, &gcnt, p*3, &gmin, p*3, &gmax));
    for (PetscInt i = 0; i < p*3; i++) { cmin[i] = PETSC_MAX_REAL; cmax[i] = PETSC_MIN_REAL; }
    for (PetscInt i = 0; i < lsz_part; i++)
      for (PetscInt s = 0; s < p; s++)
        if (Dloc[i*p + s] <= delta_s[s]) {
          ccnt[s]++;
          for (PetscInt d = 0; d < 3; d++) {
            PetscReal r = g.coords[3*i + d];
            csum[s*3 + d] += r;
            cmin[s*3 + d] = PetscMin(cmin[s*3 + d], r);
            cmax[s*3 + d] = PetscMax(cmax[s*3 + d], r);
          }
        }
    PetscCallMPI(MPI_Allreduce(csum, gsum, p*3, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(ccnt, gcnt, p,   MPIU_INT,  MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(cmin, gmin, p*3, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(cmax, gmax, p*3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PetscMalloc2(p*3, &cen, p*3, &ext));
    for (PetscInt s = 0; s < p; s++)
      for (PetscInt d = 0; d < 3; d++) {
        cen[s*3 + d] = gcnt[s] > 0 ? gsum[s*3 + d] / gcnt[s] : 0;
        ext[s*3 + d] = gmax[s*3 + d] - gmin[s*3 + d];
      }
    PetscCall(PetscFree4(csum, ccnt, cmin, cmax));
    PetscCall(PetscFree4(gsum, gcnt, gmin, gmax));
  }

  // coarse basis: smooth distance partition of unity phi_s(i) = ramp(d_s(i)) / sum_s' ramp(d_s'(i)),
  // ramp(d) = 1 - d/delta (1 at the core, 0 at the overlap edge). Bounded gradient ~1/delta keeps the
  // coarse operator well scaled. Built entirely from owned rows (no cross-rank assembly needed since
  // each node knows its distance to every subdomain). Optional enrichment columns reuse the same phi.
  ncb = nmemb * data->cb_ncomp;
  PetscCall(MatCOO_Alloc(coo, ncb));
  for (PetscInt i = 0; i < lsz_part; i++) {
    PetscReal sumr = 0;
    for (PetscInt s = 0; s < p; s++)
      if (Dloc[i*p + s] <= delta_s[s]) sumr += 1. - Dloc[i*p + s] / delta_s[s];
    if (sumr <= 0) continue; // unreachable from any core (own core has d=0, so normally impossible)
    for (PetscInt s = 0; s < p; s++) {
      if (Dloc[i*p + s] > delta_s[s]) continue;
      PetscReal phi = (1. - Dloc[i*p + s] / delta_s[s]) / sumr;
      PetscInt cb_idx = data->cb_ncomp * s;
      PetscCall(MatCOO_Push(coo, vstart+i, cb_idx, phi));                       // constant (always)
      if (cen) { // enrichment: cen/ext are allocated only when pux_dim>0 || pu_xy
        // normalized local coordinates in [-~1/2, ~1/2] (0 if the subdomain is flat along that axis)
        PetscReal nc[3];
        for (PetscInt d = 0; d < 3; d++)
          nc[d] = ext[s*3 + d] > data->eps ? (g.coords[3*i + d] - cen[s*3 + d]) / ext[s*3 + d] : 0;
        for (PetscInt d = 0; d < data->pux_dim; d++)                            // linear x, y[, z]
          PetscCall(MatCOO_Push(coo, vstart+i, cb_idx + d + 1, nc[d] * phi));
        if (data->pu_xy)                                                        // bilinear cross xy
          PetscCall(MatCOO_Push(coo, vstart+i, cb_idx + data->pux_dim + 1, nc[0] * nc[1] * phi));
      }
    }
  }

  if (data->pux_dim > 0 || data->pu_xy) PetscCall(PetscFree2(cen, ext));
  PetscCall(PetscFree(Dloc));
  PetscCall(PetscFree2(diam_s, delta_s));
  PetscCall(PetscFree(leafpart));
  PetscCall(net2as_wgraph_destroy(&g));

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
  if (!data->domain_set) PetscCall(PCSetup_Net2AS_ReadDomain(pc, comm));
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
  // The coarse basis row layout must be A's node layout (= A_local_rows / bs), NOT PETSC_DECIDE:
  // with a graph-partitioned (uneven) distribution the even split would not conform to A in the
  // MatPtAP below. (It only worked previously because round-robin produced an even split.)
  PetscInt cb_local_rows;
  PetscCall(MatGetLocalSize(A, &cb_local_rows, NULL));
  cb_local_rows /= data->bs;
  PetscCall(MatSetSizes(coarse_basis, cb_local_rows, PETSC_DECIDE, size, n_cols*data->cb_ncomp));
  PetscCall(MatSetOptionsPrefix(coarse_basis, "net2as_coarse_"));
  PetscCall(MatSetPreallocationCOO(coarse_basis, coo.nnz, coo.rows, coo.cols));
  PetscCall(MatSetValuesCOO(coarse_basis, coo.vals, INSERT_VALUES));
  PetscCall(MatCreateMAIJ(coarse_basis, data->bs, &data->cb)); // expanded by block size

  // setup coarse global data structures
  PetscCall(MatPtAP(A, data->cb, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &data->cmat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "net2as:\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  cb_type: %s\n", data->cb_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  cb_trim: %s\n", data->cb_trim ? "true" : "false"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  load_type: %s\n", data->load_type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  nocoarse: %s\n", data->nocoarse ? "true" : "false"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  p: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", data->p[0], data->p[1]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  H: %.5e\n", data->H));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  delta: %.5e\n", data->overlap_abs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  delta_rel: %.5e\n", data->overlap_frac));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  bs: %" PetscInt_FMT "\n", data->bs));
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
  // DBG probe (coarse-correction magnitude): commented out — it must NOT print mid-solve, because with
  // -ksp_norm_type unpreconditioned the first apply lands between ksp_monitor it=0 and it=1 and would
  // corrupt the YAML. Re-enable all four DBG blocks together to inspect ||A0^-1 cb^T x|| etc.
  // static int dbg_count = 0;
  // PetscReal dbg_nx = 0, dbg_nrhs = 0, dbg_ncsol = 0, dbg_ncoarse = 0, dbg_ntotal = 0;

  PetscFunctionBegin;
  if (!data->ksp) PetscCall(PCSetup_Net2AS(pc));

  PetscCall(VecScatterBegin(data->rank_sc, x, data->rank_sol, INSERT_VALUES, SCATTER_FORWARD));

  // coarse
  if (!data->nocoarse) {
    PetscCall(MatMultTranspose(data->cb, x, data->csol));
    // if (dbg_count < 1) PetscCall(VecNorm(data->csol, NORM_2, &dbg_nrhs)); // DBG: ||cb^T x||
    PetscCall(KSPSolve(data->cksp, data->csol, data->csol));
    // if (dbg_count < 1) PetscCall(VecNorm(data->csol, NORM_2, &dbg_ncsol)); // DBG: ||A0^-1 cb^T x||
    PetscCall(MatMult(data->cb, data->csol, y));
    // if (dbg_count < 1) { // DBG: ||coarse correction|| and ||x||
    //   PetscCall(VecNorm(y, NORM_2, &dbg_ncoarse));
    //   PetscCall(VecNorm(x, NORM_2, &dbg_nx));
    // }
  } else {
    // y is only otherwise initialized by the coarse MatMult above; the local correction is added
    // into it via the reverse scatter below. With no coarse correction we must zero it ourselves,
    // since PETSc does not zero the PCApply output vector on entry.
    PetscCall(VecZeroEntries(y));
  }

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

  // if (dbg_count < 1 && !data->nocoarse) { // DBG: report coarse-correction magnitude on first apply
  //   PetscCall(VecNorm(y, NORM_2, &dbg_ntotal));
  //   PetscCall(PetscPrintf(PETSC_COMM_WORLD,
  //     "net2as_dbg:\n  nx: %.5e\n  n_cb_t_x: %.5e\n  n_csol: %.5e\n  n_coarse_y: %.5e\n  n_total_y: %.5e\n  coarse_frac: %.5e\n",
  //     dbg_nx, dbg_nrhs, dbg_ncsol, dbg_ncoarse, dbg_ntotal, dbg_ntotal > 0 ? dbg_ncoarse/dbg_ntotal : 0));
  //   dbg_count++;
  // }

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
  data->overlap_frac = 0.1;
  data->overlap_abs = 0.;
  data->pux_dim = 0;
  data->pu_xy = PETSC_FALSE;
  data->cb_ncomp = 1;
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
