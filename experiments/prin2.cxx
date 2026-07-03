#include "prin2.hxx"

static PetscInt PETSC_PRIN2_ROW_LEN = 10;
static PetscLogDouble PETSC_PRIN2_TIMER = 0;
static PetscLogDouble PETSC_PRIN2_LAST_T = 0;
static PetscInt PETSC_PRIN2_STAGE = 0;
static const char* PETSC_PRIN2_STAGE_NAME = "";

PetscErrorCode PRIN2S(PetscLogStage stage) {
  PETSC_PRIN2_STAGE = stage;
  PetscCall(PetscLogStageGetName(stage, &PETSC_PRIN2_STAGE_NAME));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# %s...\n", PETSC_PRIN2_STAGE_NAME));
  PetscCall(PetscTime(&PETSC_PRIN2_TIMER));
  PetscCall(PetscLogStagePush(stage));
  return 0;
}

PetscErrorCode PRIN2SP() {
  PetscLogDouble time;
  PetscCall(PetscLogStagePop());
  PetscCall(PetscTime(&time));
  PETSC_PRIN2_LAST_T = time-PETSC_PRIN2_TIMER;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t_%s: %.5e\n",
    PETSC_PRIN2_STAGE_NAME, PETSC_PRIN2_LAST_T));
  return 0;
}

PetscErrorCode PRIN2SP(PetscLogDouble *t) {
  *t = PETSC_PRIN2_LAST_T;
  return 0;
}

PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, const PetscReal* dat, PetscInt len) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(com, msg));
  for (PetscInt i = 0; i < len; i++) {
    if (i % PETSC_PRIN2_ROW_LEN == 0)
       PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % .5e", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2i(MPI_Comm com, const char* msg, const PetscInt* dat, PetscInt len) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(com, msg));
  for (PetscInt i = 0; i < len; i++) {
    if (i % PETSC_PRIN2_ROW_LEN == 0)
      PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % 12d", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2iy(MPI_Comm comm, const char *name, PetscInt val) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: %d\n", name, val));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2fy(MPI_Comm comm, const char *name, PetscReal val) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: %.5e\n", name, val));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2iya(MPI_Comm comm, const char *name, PetscInt val, PetscInt num) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(comm, "%s: [", name));
  for (PetscInt i = 0; i < num; i++)
    PetscCall(PetscPrintf(comm, "%d, ", val));
  PetscCall(PetscPrintf(comm, "]\n"));
  PetscFunctionReturn(0);
}

PetscErrorCode PetscPrin2Options() {
  PetscBool set;

  PetscFunctionBeginUser;
  PetscOptionsBegin(PETSC_COMM_WORLD, "prin2_", "Prin2", NULL);
  PetscCall(PetscOptionsInt("-row_len", "length of displayed rows", NULL, PETSC_PRIN2_ROW_LEN, &PETSC_PRIN2_ROW_LEN, &set));
  PetscOptionsEnd();
  PetscFunctionReturn(0);
}

// must call VecRestoreSpan(x, span) after
PetscErrorCode VecGetSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p;
  PetscInt n;
  PetscCall(VecGetArray(x, &p));
  PetscCall(VecGetLocalSize(x, &n));
  span = {p, (size_t)n};
  return 0;
}

// must be called after each VecGetSpan(x, span)
PetscErrorCode VecRestoreSpan(Vec x, std::span<PetscScalar>& span) {
  PetscScalar* p = span.data();
  PetscCall(VecRestoreArray(x, &p));
  span = std::span<PetscScalar>();
  return 0;
}

PetscErrorCode KSPMonitorYAML(KSP ksp, PetscInt it, PetscReal rnorm, PetscViewerAndFormat *vf) {
  PetscViewer viewer = vf->viewer;
  KSPMonitorYAML_Ctx *ctx = (KSPMonitorYAML_Ctx*)vf->data;
  PetscLogDouble t1;

  PetscFunctionBegin;
  if (it == 0) {
    PetscCall(PetscViewerASCIIPrintf(viewer, "ksp_monitor:\n"));
  }
  PetscCall(PetscTime(&t1));
  PetscCall(PetscViewerASCIIPrintf(viewer, "  - it: %3" PetscInt_FMT "\n", it));
  PetscCall(PetscViewerASCIIPrintf(viewer, "    time: %.16e\n", (double)(t1-ctx->t0)));
  PetscCall(PetscViewerASCIIPrintf(viewer, "    rnorm: %.16e\n", (double)rnorm));
  PetscFunctionReturn(0);
}
