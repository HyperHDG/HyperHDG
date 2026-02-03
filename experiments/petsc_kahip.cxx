#include <petsc/private/matimpl.h>
#include <../src/mat/impls/adj/mpi/mpiadj.h>
#include <KaHIP/interface/kaHIP_interface.h>

struct MatPartitioning_KaHIP {
  // HACK: this must be at the same byte offset in the struct as the same field in the parmetis struct
  PetscInt cuts;
};


PetscErrorCode MatPartitioningApply_KaHIP(MatPartitioning part, IS *partition) {
  MatPartitioning_KaHIP *data = (MatPartitioning_KaHIP *)part->data;
  Mat            adj;
  PetscInt       n, *xadj, *adjncy, *adjcwgt, *parts, nparts = part->n, mode = 0, comm_size, seed = 0;
  PetscReal      imbalance = 0.03;
  PetscBool      done;
  bool suppress_output = true;

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
  kaffpa(&n, part->vertex_weights, (int*)xadj, (int*)adjcwgt, (int*)adjncy, &nparts, &imbalance, suppress_output, seed, mode, &data->cuts, parts);
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
  part->data         = (void*)ctx;
  part->ops->apply   = MatPartitioningApply_KaHIP;
  part->ops->destroy = MatPartitioningDestroy_KaHIP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

