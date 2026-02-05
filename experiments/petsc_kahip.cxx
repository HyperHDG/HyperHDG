#include <petsc/private/matimpl.h>
#include <../src/mat/impls/adj/mpi/mpiadj.h>
#include <KaHIP/parallel/parallel_src/interface/parhip_interface.h>

#define PetscArraycpyCast(dst, src, n, dsttype, srctype) \
  do { for (typeof(n) _i = 0; _i < (n); _i++) (dst)[_i] = (dsttype)((srctype*)(src))[_i]; } while (0)

struct MatPartitioning_KaHIP {
  // HACK: this must be at the same byte offset in the struct as the same field in the parmetis struct
  PetscInt cuts;
  PetscInt seed;
  PetscInt mode;
  PetscReal imbalance;
  PetscBool suppress_output;
};

PetscErrorCode MatPartitioningSetFromOptions_KaHIP(MatPartitioning part, PetscOptionItems PetscOptionsObject) {
  MatPartitioning_KaHIP *ctx = (MatPartitioning_KaHIP*)part->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject, "KaHIP Partitioning Options");
  PetscCall(PetscOptionsReal("-kahip_imbalance", "partition imbalance", NULL, ctx->imbalance, &ctx->imbalance, NULL));
  PetscCall(PetscOptionsInt("-kahip_seed", "random seed", NULL, ctx->seed, &ctx->seed, NULL));
  PetscCall(PetscOptionsBool("-kahip_suppress_output", "suppress kahip logging", NULL, ctx->suppress_output, &ctx->suppress_output, NULL));
  PetscCall(PetscOptionsInt("-kahip_mode", "kahip configuration mode", NULL, ctx->mode, &ctx->mode, NULL));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningApply_KaHIP(MatPartitioning part, IS *partition) {
  MatPartitioning_KaHIP *data = (MatPartitioning_KaHIP *)part->data;
  Mat            adj;
  PetscInt       n, m, *p_parts;
  const PetscInt *p_vtxdist, *p_xadj, *p_adjncy, *p_adjcwgt;
  PetscBool      done;

  MPI_Comm comm = PetscObjectComm((PetscObject)part);
  idxtype *vtxdist, *xadj, *adjncy, *adjcwgt, *parts, *vtxwgt;
  int seed = (int)data->seed, mode = (int)data->mode, edgecut, nparts = (int)part->n, comm_size;
  double imbalance = (double)data->imbalance;
  bool suppress = (bool)data->suppress_output;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &comm_size));
  PetscCall(MatConvert(part->adj, MATMPIADJ, MAT_INITIAL_MATRIX, &adj));
  PetscCall(MatGetOwnershipRanges(adj, &p_vtxdist)); // of size ranks+1
  PetscCall(MatGetRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, &p_xadj, &p_adjncy, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ failed");
  p_adjcwgt = ((Mat_MPIAdj*)adj->data)->values;
  m = p_xadj[n];

  // convert PetscInt to idxtype
  PetscCheck(sizeof(PetscInt) <= sizeof(idxtype), PETSC_COMM_SELF, PETSC_ERR_PLIB,
    "ParHIP: sizeof(PetscInt) == %lu must be at most sizeof(idxtype) == %lu", sizeof(PetscInt), sizeof(idxtype));
  PetscCall(PetscMalloc7(n+1, &xadj, m, &adjncy, m, &adjcwgt, n, &parts, n, &p_parts, n, &vtxwgt, comm_size+1, &vtxdist));
  PetscArraycpyCast(xadj, p_xadj, n+1, idxtype, PetscInt);
  PetscArraycpyCast(adjncy, p_adjncy, m, idxtype, PetscInt);
  (void)p_adjcwgt;
  PetscArraycpyCast(vtxdist, p_vtxdist, comm_size+1, idxtype, PetscInt);

  // perform partitioning
  ParHIPPartitionKWay(vtxdist, xadj, adjncy, NULL, NULL,
    &nparts, &imbalance, suppress, seed, mode, &edgecut, parts, &comm);

  // NOTE: narrowing cast is ok as the number of partitions already is a PetscInt
  data->cuts = edgecut;
  PetscArraycpyCast(p_parts, parts, n, PetscInt, idxtype);
  PetscCall(ISCreateGeneral(comm, n, p_parts, PETSC_COPY_VALUES, partition));

  // NOTE: zeros n
  PetscCall(MatRestoreRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, &p_xadj, &p_adjncy, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ failed");
  PetscCall(PetscFree7(xadj, adjncy, adjcwgt, parts, p_parts, vtxwgt, vtxdist));
  PetscCall(MatDestroy(&adj));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningDestroy_KaHIP(MatPartitioning part) {
  PetscFunctionBegin;
  PetscCall(PetscFree(part->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningCreate_KaHIP(MatPartitioning part) {
  MatPartitioning_KaHIP *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ctx->suppress_output = PETSC_TRUE;
  ctx->seed = 0;
  ctx->imbalance = 0.03;
  ctx->mode = ULTRAFASTMESH; // 0
  part->data         = (void*)ctx;
  part->ops->setfromoptions = MatPartitioningSetFromOptions_KaHIP;
  part->ops->apply   = MatPartitioningApply_KaHIP;
  part->ops->destroy = MatPartitioningDestroy_KaHIP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

