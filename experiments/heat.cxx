#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_parab_ldgh.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include "parameters.hxx"
#include "../reproducibles_python/parameters/diffusion.hxx"

static const char help[] = "experiments regarding the heat equation\n";

static const char help_msg[] = "experiments regarding timoshenko networks\n";
// static PetscInt PETSC_PRIN2_ROW_LEN = 10;
static PetscLogDouble PETSC_PRIN2_TIMER = 0;
static PetscInt PETSC_PRIN2_STAGE = 0;
static const char* PETSC_PRIN2_STAGE_NAME = "";

#define PRIN2IY(VAR)  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2FY(VAR)  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2SY(VAR)  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", #VAR, VAR))
#define PRIN2S(STAGE) do { PETSC_PRIN2_STAGE = STAGE; PetscCall(PetscLogStageGetName(STAGE, &PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# %s...\n", PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscTime(&PETSC_PRIN2_TIMER)); PetscCall(PetscLogStagePush(STAGE)); } while(0)
#define PRIN2SP()     do { PetscLogDouble time; PetscCall(PetscLogStagePop()); PetscCall(PetscTime(&time)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t_%s: %.5e\n", PETSC_PRIN2_STAGE_NAME, (time-PETSC_PRIN2_TIMER))); } while(0)

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

PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, PetscReal* dat, PetscInt len) {
  PetscCall(PetscPrintf(com, msg));
  const PetscInt row_len = 10;
  for (PetscInt i = 0; i < len; i++) {
    if (i % row_len == 0)
       PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % .5e", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
  return 0;
}

PetscErrorCode PetscPrin2i(MPI_Comm com, const char* msg, PetscInt* dat, PetscInt len) {
  PetscCall(PetscPrintf(com, msg));
  const PetscInt row_len = 10;
  for (PetscInt i = 0; i < len; i++) {
    if (i % row_len == 0)
      PetscCall(PetscPrintf(com, "\n"));
    PetscCall(PetscPrintf(com, "  % 12d", dat[i]));
  }
  PetscCall(PetscPrintf(com, "\n"));
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

int main(int argc, char **argv) {
    constexpr int space_dim = 1;
    constexpr int poly_deg = 3;
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionParab<space_dim,poly_deg,2*poly_deg,TestParametersSinParab,PetscReal>;
    using HDG = GlobalLoop::Parabolic<Top,Geo,NDes,LSol>;

    PetscReal tau = 1; // HDG penalty
    PetscReal theta = .5; // one-step theta method
    PetscInt it = 2;
    PetscInt timesteps = 100;
    PetscReal end_time = 1;
    PetscReal dt;
    PetscReal rtol = 1e-13;

    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "heat";

    PetscLogStage s_as, s_ts, s_rf, s_hdg;

    PetscBool is_set;
    PetscInt N;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, e_abs = 0, e_rel = 0;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    std::vector<PetscInt> itemp;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "# initialization...\n"));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-theta", &theta, &is_set));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-i", &it, &is_set));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-ts", &timesteps, &is_set));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &end_time, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-o", output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-od", output_directory, PATH_MAX, &is_set));
    PetscCall(PetscLogStageRegister("assembly", &s_as));
    PetscCall(PetscLogStageRegister("timestepping", &s_ts));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));
    PetscCall(PetscLogStageRegister("hdg_init", &s_hdg));

    timesteps = 1<<timesteps;
    dt = end_time / timesteps;

    PRIN2IY(space_dim);
    PRIN2IY(poly_deg);
    PRIN2FY(tau);
    PRIN2FY(theta);
    PRIN2IY(it);
    PRIN2IY(tau);
    PRIN2IY(timesteps);
    PRIN2FY(dt);
    PRIN2FY(end_time);

    PRIN2S(s_hdg);
    HDG hdg((1 << it) * space_dim, {tau, theta, dt});
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "true");
    hdg.plot_option("scale", "0.95");

    zero_v = hdg.zero_vector();
    temp = hdg.make_initial(zero_v);
    N = temp.size();
    hdg.plot_solution(temp, 0.); // needs petsc
    PRIN2SP();

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));

    PRIN2S(s_as);
    PetscCall(MatCreateFromOptions(PETSC_COMM_WORLD, "t2f_", 1, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    mat_coo = hdg.trace_to_flux_mat(0.);
    PetscCall(MatSetPreallocationCOO(mat, mat_coo.row_vec.size(), (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, (PetscReal*)mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, /* keep = */ PETSC_FALSE));

    PRIN2SP();

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPSetFromOptions(ksp));

    Vec errors;
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, timesteps, timesteps, &errors));

    PRIN2S(s_ts);
    for (PetscInt i = 0; i < timesteps; i++) {
        std::span<PetscReal> rhs_span;
        std::span<PetscReal> sol_span;
        PetscCall(VecGetSpan(rhs, rhs_span));
        PetscCall(VecGetSpan(sol, sol_span));

        PetscLogStagePush(s_rf);
        hdg.residual_flux2(std::span{zero_v}, rhs_span, (i+1)*dt);
        PetscLogStagePop();
        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, sol));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        hdg.set_data(sol_span, (i+1)*dt);
        hdg.plot_solution(sol_span, (i+1)*dt);

        temp2 = hdg.errors(temp, (i+1)*dt);
        temp3 = hdg.norms(temp, (i+1)*dt);
        e_abs = PetscMax(temp2[0], e_abs);
        e_rel = PetscMax(temp2[0] / temp3[0], e_rel);
        PetscCall(VecSetValue(errors, i, temp2[0]/temp3[0], INSERT_VALUES));

        PetscCall(VecRestoreSpan(rhs, rhs_span));
        PetscCall(VecRestoreSpan(sol, sol_span));
    }
    PRIN2SP();

    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));

    avg_iterations = ((PetscReal)iterations) / timesteps;

    PRIN2FY(e_abs);
    PRIN2FY(e_rel);
    PRIN2IY(iterations);
    PRIN2FY(avg_iterations);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "output: '%s/%s.*.vtu'\n", output_directory, output_filename));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}
