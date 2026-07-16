#include "petsc_util.hxx"
#include "petsc.h"

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

PetscErrorCode MatPrintSymmetry(const char* msg, Mat mat) {
  Mat AT, D;
  PetscReal nrm, nrm_a;

  PetscFunctionBeginUser;
  MatTranspose(mat, MAT_INITIAL_MATRIX, &AT);
  MatDuplicate(mat, MAT_COPY_VALUES, &D);
  MatAXPY(D, -1.0, AT, DIFFERENT_NONZERO_PATTERN);
  MatNorm(D, NORM_FROBENIUS, &nrm);
  MatNorm(mat, NORM_FROBENIUS, &nrm_a);
  PetscPrintf(PETSC_COMM_WORLD, "%s: %g\n", msg, (double)(nrm/nrm_a));
  MatDestroy(&AT); MatDestroy(&D);
  PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsLeftYAML(PetscOptions options) {
    PetscInt unused;
    char **names;
    char **values;

    PetscCall(PetscOptionsLeftGet(NULL, &unused, &names, &values));
    if (unused == 0) goto end;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# WARNING! There are options you set that were not used!\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "options_left:\n"));
    for (PetscInt i = 0; i < unused; i++)
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  - name: \"%s\"\n    value: \"%s\"\n", names[i], values[i]));
end:
    PetscCall(PetscOptionsLeftRestore(NULL, &unused, &names, &values));
    PetscCall(PetscOptionsSetValue(NULL, "-options_left", "0"));
    return 0;
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

PetscErrorCode KSPMonitorYAML_Setup(KSP ksp, Vec rhs, void *ctx_) {
  KSPMonitorYAML_Ctx *ctx = (KSPMonitorYAML_Ctx*)ctx_;
  PetscBool is_set;
  PetscReal rtol_user, t0, t_ref, rtol_ref = 1e-14;
  PetscInt ref_its;
  const char* ref_creason;
  Mat mat;
  Vec sol_ref;

  ctx->enorm = false;

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "KSPMonitorYAML Options", NULL);
  PetscCall(PetscOptionsBool("-ksp_monitor_yaml_enorm", "compute energy norm error", NULL, ctx->enorm, &ctx->enorm, &is_set));
  PetscCall(PetscOptionsReal("-ksp_monitor_yaml_enorm_rtol_ref", "residual tolerance for reference solution used in energy error", NULL, rtol_ref, &rtol_ref, &is_set));
  PetscOptionsEnd();

  if (!ctx->enorm) goto end;

  ctx->quiet = PETSC_TRUE;
  PetscCall(KSPGetOperators(ksp, &mat, &mat));
  PetscCall(KSPGetTolerances(ksp, &rtol_user, NULL, NULL, NULL));
  PetscCall(MatCreateVecs(mat, &sol_ref, NULL));
  PetscCall(KSPSetTolerances(ksp, rtol_ref, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
  PetscCall(PetscTime(&t0));
  PetscCall(KSPSolve(ksp, rhs, sol_ref));
  PetscCall(PetscTime(&t_ref));
  t_ref -= t0;
  // exclude the reference solve from the main solve's iteration timestamps
  ctx->t0 += t_ref;
  PetscCall(KSPGetConvergedReasonString(ksp, &ref_creason));
  PetscCall(KSPGetIterationNumber(ksp, &ref_its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "ref_its: %" PetscInt_FMT "\n", ref_its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "ref_creason: %s\n", ref_creason));
  PetscCall(KSPSetTolerances(ksp, rtol_user, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
  ctx->mat = mat;
  ctx->u_ref = sol_ref;
  PetscCall(VecDuplicate(sol_ref, &ctx->e));
  PetscCall(VecDuplicate(sol_ref, &ctx->Ke));
  ctx->quiet = PETSC_FALSE;

  // installed only now so the reference solve above still uses the default residual test
  PetscCall(KSPSetConvergenceTest(ksp, KSPConvergedEnorm, ctx, NULL));

end:
  return 0;
}

// |u_ref - u^(it)|_K for the current iterate; caches the value in ctx so the convergence test can
// reuse the monitor's computation within the same iteration
static PetscErrorCode KSPEnormCompute(KSP ksp, PetscInt it, KSPMonitorYAML_Ctx *ctx, PetscReal *enorm) {
  Vec sol;
  PetscScalar dot;

  PetscFunctionBeginUser;
  PetscCall(KSPBuildSolution(ksp, NULL, &sol));
  PetscCall(VecWAXPY(ctx->e, -1.0, sol, ctx->u_ref));
  PetscCall(MatMult(ctx->mat, ctx->e, ctx->Ke));
  PetscCall(VecDot(ctx->e, ctx->Ke, &dot));
  *enorm = PetscSqrtReal(PetscMax(PetscRealPart(dot), 0.));
  ctx->enorm_it = it;
  ctx->enorm_val = *enorm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Energy-error convergence test: stop when |u_ref - u^(it)|_K <= rtol * |u_ref - u^(0)|_K
// (KSPConvergedDefault semantics with the energy error in place of the residual norm). Needs the
// reference solution from KSPMonitorYAML_Setup; enabled with -ksp_converged_enorm.
PetscErrorCode KSPConvergedEnorm(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason *reason, void *ctx_) {
  KSPMonitorYAML_Ctx *ctx = (KSPMonitorYAML_Ctx*)ctx_;
  PetscReal enorm, rtol, abstol, dtol;

  PetscFunctionBeginUser;
  (void)rnorm;
  *reason = KSP_CONVERGED_ITERATING;
  if (ctx->enorm_it == it) enorm = ctx->enorm_val;
  else PetscCall(KSPEnormCompute(ksp, it, ctx, &enorm));
  // consume the cache so a later solve restarting at the same it cannot reuse a stale value
  ctx->enorm_it = -1;
  if (it == 0) ctx->enorm0 = enorm;
  PetscCall(KSPGetTolerances(ksp, &rtol, &abstol, &dtol, NULL));
  if (PetscIsInfOrNanReal(enorm)) *reason = KSP_DIVERGED_NANORINF;
  else if (enorm <= rtol * ctx->enorm0) *reason = KSP_CONVERGED_RTOL;
  else if (enorm <= abstol) *reason = KSP_CONVERGED_ATOL;
  else if (enorm > dtol * ctx->enorm0) *reason = KSP_DIVERGED_DTOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Per-iteration energy-norm error |u_ref - u^(it)|_K against a converged reference solution
// (gortz.pdf fig. 8).
PetscErrorCode KSPMonitorYAML(KSP ksp, PetscInt it, PetscReal rnorm, PetscViewerAndFormat *vf) {
  // KSPMonitorSetFromOptions stashes the user ctx in vf->data
  KSPMonitorYAML_Ctx *ctx = (KSPMonitorYAML_Ctx*)vf->data;
  PetscReal enorm;
  PetscLogDouble t1;

  PetscFunctionBeginUser;
  if (ctx->quiet) PetscFunctionReturn(PETSC_SUCCESS);
  if (it == 0) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "ksp_monitor:\n"));
  if (ctx->enorm) PetscCall(KSPEnormCompute(ksp, it, ctx, &enorm));
  PetscCall(PetscTime(&t1));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  - it: %3" PetscInt_FMT "\n", it));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    time: %.16e\n", (double)(t1 - ctx->t0)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    rnorm: %.16e\n", (double)rnorm));
  if (ctx->enorm)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    enorm: %.16e\n", (double)enorm/ctx->enorm0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

