#pragma once

#include <petsc.h>

struct PC_Net2AS {
  // path to domain file
  char domain[PATH_MAX];
  // flat coordinate array in row-major ordering, x0,y0,z0,x1,...
  Vec points;
  // bounding box of points
  PetscReal min[3], max[3];
  // number of subdomains in [x,y]
  PetscInt p[2];
  // number of local data structures
  PetscInt sz;
  // block size (number of dofs per node)
  PetscInt bs;

  // coarse basis representation of the overlapping subdomains,
  // expanded by block size
  Mat  cb;

  // local datastructures
  // layout: i=0 -> coarse, 1 <= i < sz -> local, corresponding to subdomains
  //   ignore is[0]
  Mat* mat;
  KSP* ksp;
  IS* is;
  Vec* sol;
  VecScatter* sc;
};

PetscErrorCode PCDestroy_Net2AS(PC pc);
PetscErrorCode PCSetFromOptions_Net2AS(PC pc, PetscOptionItems PetscOptionsObject);
PetscErrorCode PCSetup_Net2AS_ReadDomain(PC pc, MPI_Comm comm);
PetscErrorCode PCSetup_Net2AS_SetupKSP(PC pc, MPI_Comm comm, PetscInt i);
PetscErrorCode PCSetup_Net2AS(PC pc);
PetscErrorCode PCApply_Net2AS(PC pc, Vec x, Vec y);
PetscErrorCode PCView_Net2AS(PC pc, PetscViewer viewer);
PetscErrorCode PCCreate_Net2AS(PC pc);
