#include <petsc.h>
#include <span>

#define PRIN2IY(VAR)  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2FY(VAR)  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2SY(VAR)  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", #VAR, VAR))

PetscErrorCode PRIN2S(PetscLogStage stage);
PetscErrorCode PRIN2SP();
PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, const PetscReal* dat, PetscInt len);
PetscErrorCode PetscPrin2i(MPI_Comm com, const char* msg, const PetscInt* dat, PetscInt len);
PetscErrorCode PetscPrin2iy(MPI_Comm comm, const char *name, PetscInt val);
PetscErrorCode PetscPrin2fy(MPI_Comm comm, const char *name, PetscReal val);
PetscErrorCode PetscPrin2iya(MPI_Comm comm, const char *name, PetscInt val, PetscInt num);
PetscErrorCode PetscPrin2Options();
// must call VecRestoreSpan(x, span) after
PetscErrorCode VecGetSpan(Vec x, std::span<PetscScalar>& span);
// must be called after each VecGetSpan(x, span)
PetscErrorCode VecRestoreSpan(Vec x, std::span<PetscScalar>& span);
