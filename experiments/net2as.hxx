#pragma once

#include <petsc.h>

PetscErrorCode PCDestroy_Net2AS(PC pc);
PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject);
PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm);
PetscErrorCode PCSetup_Net2AS_SetupKSP(PC pc, MPI_Comm comm, PetscInt i);
PetscErrorCode PCSetup_Net2AS(PC pc);
PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y);
PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer);
PetscErrorCode PCCreate_Net2AS(PC pc);
PetscErrorCode PCNet2ASGetCB(PC pc, Mat *cb);
