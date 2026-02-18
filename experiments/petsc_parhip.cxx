#include <petsc/private/matimpl.h>
#include <KaHIP/parallel/parallel_src/interface/parhip_interface.h>

#define PetscArraycpyCast(dst, src, n, dsttype, srctype) \
  do { for (typeof(n) _i = 0; _i < (n); _i++) (dst)[_i] = (dsttype)((srctype*)(src))[_i]; } while (0)

struct MatPartitioning_ParHIP {
  // HACK: this must be at the same byte offset in the struct as the same field in the parmetis struct
  PetscInt cuts;
  PetscInt seed;
  PetscInt mode;
  PetscReal imbalance;
  PetscBool suppress_output;
};

PetscErrorCode MatPartitioningSetFromOptions_ParHIP(MatPartitioning part, PetscOptionItems PetscOptionsObject) {
  MatPartitioning_ParHIP *ctx = (MatPartitioning_ParHIP*)part->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject, "ParHIP Partitioning Options");
  PetscCall(PetscOptionsReal("-parhip_imbalance", "partition imbalance", NULL, ctx->imbalance, &ctx->imbalance, NULL));
  PetscCall(PetscOptionsInt("-parhip_seed", "random seed", NULL, ctx->seed, &ctx->seed, NULL));
  PetscCall(PetscOptionsBool("-parhip_suppress_output", "suppress parhip logging", NULL, ctx->suppress_output, &ctx->suppress_output, NULL));
  PetscCall(PetscOptionsInt("-parhip_mode", "parhip configuration mode", NULL, ctx->mode, &ctx->mode, NULL));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningApply_ParHIP(MatPartitioning part, IS *partition) {
  MatPartitioning_ParHIP *data = (MatPartitioning_ParHIP *)part->data;
  Mat            adj;
  PetscInt       n, m, *p_parts;
  const PetscInt *p_vtxdist, *p_xadj, *p_adjncy;
  PetscBool      done;

  MPI_Comm comm = PetscObjectComm((PetscObject)part);
  idxtype *vtxdist, *xadj, *adjncy, *adjcwgt = NULL, *parts, *vtxwgt;
  int seed = (int)data->seed, mode = (int)data->mode, edgecut, nparts = (int)part->n, comm_size;
  double imbalance = (double)data->imbalance;
  bool suppress = (bool)data->suppress_output;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &comm_size));
  PetscCall(MatConvert(part->adj, MATMPIADJ, MAT_INITIAL_MATRIX, &adj));
  PetscCall(MatGetOwnershipRanges(adj, &p_vtxdist)); // of size ranks+1
  PetscCall(MatGetRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, &p_xadj, &p_adjncy, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ failed");
  m = p_xadj[n];

  // convert PetscInt to idxtype
  PetscCheck(sizeof(PetscInt) <= sizeof(idxtype), PETSC_COMM_SELF, PETSC_ERR_PLIB,
    "ParHIP: sizeof(PetscInt) == %lu must be at most sizeof(idxtype) == %lu", sizeof(PetscInt), sizeof(idxtype));
  PetscCall(PetscMalloc6(n+1, &xadj, m, &adjncy, n, &parts, n, &p_parts, n, &vtxwgt, comm_size+1, &vtxdist));
  PetscArraycpyCast(xadj, p_xadj, n+1, idxtype, PetscInt);
  PetscArraycpyCast(adjncy, p_adjncy, m, idxtype, PetscInt);
  PetscArraycpyCast(vtxwgt, part->vertex_weights, n, idxtype, PetscInt);
  PetscArraycpyCast(vtxdist, p_vtxdist, comm_size+1, idxtype, PetscInt);

  // perform partitioning
  ParHIPPartitionKWay(vtxdist, xadj, adjncy, vtxwgt, adjcwgt,
    &nparts, &imbalance, suppress, seed, mode, &edgecut, parts, &comm);

  // NOTE: narrowing cast is ok as the number of partitions already is a PetscInt
  data->cuts = edgecut;
  PetscArraycpyCast(p_parts, parts, n, PetscInt, idxtype);
  PetscCall(ISCreateGeneral(comm, n, p_parts, PETSC_COPY_VALUES, partition));

  // NOTE: zeros n
  PetscCall(MatRestoreRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, &p_xadj, &p_adjncy, &done));
  PetscCheck(done, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "MatGetRowIJ failed");
  PetscCall(PetscFree6(xadj, adjncy, parts, p_parts, vtxwgt, vtxdist));
  PetscCall(MatDestroy(&adj));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningDestroy_ParHIP(MatPartitioning part) {
  PetscFunctionBegin;
  PetscCall(PetscFree(part->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatPartitioningCreate_ParHIP(MatPartitioning part) {
  MatPartitioning_ParHIP *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ctx->suppress_output = PETSC_TRUE;
  ctx->seed = 0;
  ctx->imbalance = 0.03;
  ctx->mode = ULTRAFASTMESH; // 0
  part->data         = (void*)ctx;
  part->ops->setfromoptions = MatPartitioningSetFromOptions_ParHIP;
  part->ops->apply   = MatPartitioningApply_ParHIP;
  part->ops->destroy = MatPartitioningDestroy_ParHIP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

