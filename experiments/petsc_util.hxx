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
// print |A - A^T|_F / |A|_F
PetscErrorCode MatPrintSymmetry(const char* msg, Mat mat);
// print unused options as a YAML list, then suppress PETSc's own report
PetscErrorCode PetscOptionsLeftYAML(PetscOptions options);
// must call VecRestoreSpan(x, span) after
PetscErrorCode VecGetSpan(Vec x, std::span<PetscScalar>& span);
// must be called after each VecGetSpan(x, span)
PetscErrorCode VecRestoreSpan(Vec x, std::span<PetscScalar>& span);

struct KSPMonitorYAML_Ctx {
  // monitor fires on every KSPSolve once installed; quiet suppresses it during the
  // reference solve in KSPMonitorYAML_Setup
  PetscBool quiet = PETSC_FALSE;
  PetscBool enorm = PETSC_FALSE;
  Mat mat = NULL;
  Vec u_ref = NULL, e = NULL, Ke = NULL;
  // initial timepoint relative to which the iteration times are measured
  PetscLogDouble t0 = 0;
  // |u_ref - u^(0)|_K, fixed by KSPConvergedEnorm at it 0 of the main solve
  PetscReal enorm0 = 0;
  // energy error the monitor computed this iteration; the convergence test (which runs after the
  // monitor within an iteration) reuses it to avoid a second MatMult (-1 = no valid entry)
  PetscInt enorm_it = -1;
  PetscReal enorm_val = 0;
};
PetscErrorCode KSPMonitorYAML(KSP ksp, PetscInt it, PetscReal rnorm, PetscViewerAndFormat *vf);
PetscErrorCode KSPConvergedEnorm(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason *reason, void *ctx);
PetscErrorCode KSPMonitorYAML_Setup(KSP ksp, Vec rhs, void *ctx);
