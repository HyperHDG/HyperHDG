#pragma once

#include <petsc.h>

PetscErrorCode MatPartitioningApply_KaHIP(MatPartitioning part, IS *partition);
PetscErrorCode MatPartitioningDestroy_KaHIP(MatPartitioning part);
PetscErrorCode MatPartitioningCreate_KaHIP(MatPartitioning part);

