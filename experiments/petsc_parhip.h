#pragma once

#include <petsc.h>

PetscErrorCode MatPartitioningApply_ParHIP(MatPartitioning part, IS *partition);
PetscErrorCode MatPartitioningDestroy_ParHIP(MatPartitioning part);
PetscErrorCode MatPartitioningCreate_ParHIP(MatPartitioning part);

