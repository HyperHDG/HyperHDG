#include <petsc/private/matimpl.h>
#include <../src/mat/impls/adj/mpi/mpiadj.h>
#include <KaHIP/interface/kaHIP_interface.h>

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
  PetscInt       n, *xadj, *adjncy, *adjcwgt, *parts, nparts = part->n, comm_size;
  PetscBool      done;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)part), &comm_size));
  PetscCheck(comm_size == 1, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "kahip partitioner expects a single process");
  PetscCheck(sizeof(PetscInt) == sizeof(int), PETSC_COMM_WORLD, PETSC_ERR_PLIB, "kahip partitioner expects 32-bit PetscInt");
  PetscCheck(sizeof(PetscReal) == sizeof(double), PETSC_COMM_WORLD, PETSC_ERR_PLIB, "kahip partitioner expects 64-bit PetscReal");

  PetscCall(MatConvert(part->adj, MATMPIADJ, MAT_INITIAL_MATRIX, &adj));
  adjcwgt = ((Mat_MPIAdj*)adj->data)->values;

  PetscCall(MatGetLocalSize(adj, &n, NULL));
  PetscCall(PetscMalloc1(n, &parts));

  PetscCall(MatGetRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, (const PetscInt**)&xadj, (const PetscInt**)&adjncy, &done));
  kaffpa(&n, part->vertex_weights, (int*)xadj, (int*)adjcwgt, (int*)adjncy, &nparts, &data->imbalance, data->suppress_output, data->seed, data->mode, &data->cuts, parts);
  PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)part), n, parts, PETSC_OWN_POINTER, partition));
  // NOTE: zeros n
  PetscCall(MatRestoreRowIJ(adj, 0, PETSC_FALSE, PETSC_FALSE, &n, (const PetscInt**)&xadj, (const PetscInt**)&adjncy, &done));
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
  ctx->imbalance = 0.3;
  ctx->mode = FAST;
  part->data         = (void*)ctx;
  part->ops->setfromoptions = MatPartitioningSetFromOptions_KaHIP;
  part->ops->apply   = MatPartitioningApply_KaHIP;
  part->ops->destroy = MatPartitioningDestroy_KaHIP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

