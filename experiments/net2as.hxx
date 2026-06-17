#pragma once

#include <petsc.h>

PetscErrorCode PCDestroy_Net2AS(PC pc);
PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject);
PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm);
PetscErrorCode PCNet2ASSetDomain(PC pc,
                                 PetscInt n_owned_nodes,
                                 PetscInt n_global_nodes,
                                 PetscInt sdim,
                                 const PetscReal* coords,
                                 PetscInt n_owned_edges,
                                 const PetscInt* edges_global);
PetscErrorCode PCSetup_Net2AS_SetupKSP(PC pc, MPI_Comm comm, PetscInt i);
PetscErrorCode PCSetup_Net2AS(PC pc);
PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y);
PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer);
PetscErrorCode PCCreate_Net2AS(PC pc);
PetscErrorCode PCNet2ASGetCB(PC pc, Mat* cb);
