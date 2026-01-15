#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_wave1_ldgh.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include "parameters.hxx"
#include <map>

static const char help_msg[] = "experiments regarding the wave equation\n";
static PetscInt PETSC_PRIN2_ROW_LEN = 10;
static PetscLogDouble PETSC_PRIN2_TIMER = 0;
static PetscInt PETSC_PRIN2_STAGE = 0;
static const char* PETSC_PRIN2_STAGE_NAME = "";

#define PRIN2IY(VAR)  PetscCall(PetscPrin2iy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2FY(VAR)  PetscCall(PetscPrin2fy(PETSC_COMM_WORLD, #VAR, VAR))
#define PRIN2SY(VAR)  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", #VAR, VAR))
#define PRIN2S(STAGE) do { PETSC_PRIN2_STAGE = STAGE; PetscCall(PetscLogStageGetName(STAGE, &PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "# %s...\n", PETSC_PRIN2_STAGE_NAME)); PetscCall(PetscTime(&PETSC_PRIN2_TIMER)); PetscCall(PetscLogStagePush(STAGE)); } while(0)
#define PRIN2SP()     do { PetscLogDouble time; PetscCall(PetscLogStagePop()); PetscCall(PetscTime(&time)); PetscCall(PetscPrintf(PETSC_COMM_WORLD, "t_%s: %.5e\n", PETSC_PRIN2_STAGE_NAME, (time-PETSC_PRIN2_TIMER))); } while(0)

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



int main(int argc, char **argv) {
    constexpr int space_dim = 1;
    constexpr int poly_deg = 3;
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionWave1<space_dim,poly_deg,2*poly_deg,TestWave1,PetscReal>;
    using HDG = GlobalLoop::Hyperbolic<Top,Geo,NDes,LSol>;

    PetscBool help = false, is_set;
    PetscInt nx = 2, nt = 1;
    PetscInt N;            // global system size
    PetscReal tau = 1;     // HDG penalty
    PetscReal theta = .5;  // one-step theta method
    PetscReal T = 1, dt = 0, rtol = 1e-10, h = 0, e_abs = 0, e_rel = 0;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_iterations = 0, rnorm;
    const char* creason = NULL;
    PetscBool plot = true;
    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "wave";
    char plot_scale[PATH_MAX] = "0.95";

    PetscLogStage s_as, s_ts, s_rf;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    std::vector<PetscInt> itemp;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol, errors;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "HDG Wave Equation Options", NULL);
    PetscCall(PetscOptionsReal("-theta", "time-step averaging weight, 0 < theta <= 0.5, use theta=0.25 for CN", NULL, theta, &theta, &is_set));
    PetscCall(PetscOptionsReal("-tau", "hdg penalty parameter, recommended: tau ~ h^s for s in {-1,0,1}", NULL, tau, &tau, &is_set));
    PetscCall(PetscOptionsInt("-nx", "number of elements divide the domain into", NULL, nx, &nx, &is_set));
    PetscCall(PetscOptionsInt("-nt", "number of timesteps", NULL, nt, &nt, &is_set));
    PetscCall(PetscOptionsReal("-T", "end time", NULL, T, &T, &is_set));
    PetscCall(PetscOptionsString("-o", "output filename", NULL, output_filename, output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsString("-od", "output directory", NULL, output_directory, output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsBool("-plot", "plot solution", NULL, plot, &plot, &is_set));
    PetscCall(PetscOptionsString("-plot_scale", "subdomain scale factor for plotting", NULL, plot_scale, plot_scale, PATH_MAX, &is_set));
    PetscOptionsEnd();

    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
    if (help) {
      PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
      PetscFinalize();
      return 0;
    }

    PetscCall(PetscLogStageRegister("Assembly", &s_as));
    PetscCall(PetscLogStageRegister("Timestepping", &s_ts));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));

    dt = T / nt;
    h = 1. / nx;

    HDG hdg(nx, {tau, theta, dt});
    PRIN2IY(space_dim);
    PRIN2IY(poly_deg);
    PRIN2FY(tau);
    PRIN2FY(theta);
    PRIN2IY(nt);
    PRIN2FY(dt);
    PRIN2FY(T);
    PRIN2FY(h);
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "true");
    hdg.plot_option("scale", plot_scale);

    zero_v = hdg.zero_vector();
    N = zero_v.size();
    temp = hdg.make_initial(zero_v);
    if (plot)
      hdg.plot_solution(temp, 0.);

    temp2 = hdg.errors(temp, 0);
    temp3 = hdg.norms(temp, 0);
    e_abs = PetscMax(temp2[0], e_abs);
    e_rel = PetscMax(temp2[0] / temp3[0], e_rel);

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));
    PetscCall(VecCreateFromOptions(PETSC_COMM_SELF, "err_", 1, nt, nt, &errors));
    PetscCall(VecSetValue(errors, 0, temp2[0]/temp3[0], INSERT_VALUES));

    PetscCall(PetscPrintf(PETSC_COMM_SELF, "e_abs0: %.5e\n", e_abs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "e_rel0: %.5e\n", e_rel));

    PRIN2S(s_as);
    mat_coo = hdg.trace_to_flux_mat(0.);
    PetscCall(MatCreateFromOptions(PETSC_COMM_SELF, NULL, 1, PETSC_DECIDE, PETSC_DECIDE, N, N, &mat));
    PetscCall(MatSetPreallocationCOO(mat, mat_coo.value_vec.size(), (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_coo.value_vec.data(), INSERT_VALUES));
    PetscCall(MatEliminateZeros(mat, PETSC_TRUE));
    PRIN2SP();

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetTolerances(ksp, rtol, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    PetscCall(KSPSetFromOptions(ksp));

    PRIN2S(s_ts);
    for (PetscInt i = 1; i <= nt; i++) {
        PetscReal ti = i*dt;

        std::span<PetscReal> span;
        PetscCall(VecGetSpan(rhs, span));
        PetscLogStagePush(s_rf);
        hdg.residual_flux2(std::span{zero_v}, span, ti);
        PetscLogStagePop();

        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, rhs));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        hdg.set_data(span, ti);
        if (plot)
          hdg.plot_solution(span, ti);
        PetscCall(VecRestoreSpan(rhs, span));

        temp2 = hdg.errors(temp, ti);
        temp3 = hdg.norms(temp, ti);
        e_abs = PetscMax(temp2[0], e_abs);
        e_rel = PetscMax(temp2[0] / temp3[0], e_rel);
        PetscCall(VecSetValue(errors, i, temp2[0]/temp3[0], INSERT_VALUES));
    }
    PRIN2SP();

    PetscCall(KSPGetConvergedReasonString(ksp, &creason));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));

    PetscCall(VecAssemblyBegin(errors));
    PetscCall(VecAssemblyEnd(errors));

    avg_iterations = ((PetscReal)iterations) / nt;

    PRIN2FY(e_abs);
    PRIN2FY(e_rel);
    PRIN2IY(iterations);
    PRIN2FY(rnorm);
    PRIN2SY(creason);
    PRIN2FY(avg_iterations);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "output: %s/%s.*.vtu\n", output_directory, output_filename));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}
