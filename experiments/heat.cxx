#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/local_solver/diffusion_parab_ldgh.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include "parameters.hxx"

static const char help[] = "experiments regarding the heat equation\n";

PetscErrorCode PetscPrin2f(MPI_Comm com, const char* msg, PetscReal* dat, PetscInt len) {
  PetscCall(PetscPrintf(com, msg));
  const PetscInt row_len = 10;
  for (PetscInt i = 0; i < len; i++) {
    if (i % 2 == 0)
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

int main(int argc, char **argv) {
    int errcode = 0;

    constexpr int space_dim = 2;
    constexpr int poly_deg = 3;
    using Top = Topology::Cubic<space_dim,space_dim>;
    using Geo = Geometry::UnitCube<space_dim,space_dim,PetscReal>;
    using NDes = NodeDescriptor::Cubic<space_dim,space_dim>;
    using LSol = LocalSolver::DiffusionParab<space_dim,poly_deg,2*poly_deg,TestHeat,PetscReal>;
    using HDG = GlobalLoop::Parabolic<Top,Geo,NDes,LSol>;

    PetscReal tau = 1; // HDG penalty
    PetscReal theta = .5; // one-step theta method
    PetscInt iteration = 2;
    PetscInt timesteps = 100;
    PetscReal end_time = 1;
    PetscReal dt;

    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "heat";

    PetscLogStage stage_assembly, stage_timestep;

    PetscBool is_set;
    PetscInt N;
    PetscReal err, sol_norm;
    PetscInt iterations = 0, its = 0;
    PetscReal avg_it = 0;

    std::vector<PetscReal> temp, temp2, temp3;
    std::vector<PetscInt> itemp;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "setup...\n"));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-theta", &theta, &is_set));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-i", &iteration, &is_set));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-ts", &timesteps, &is_set));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &end_time, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-o", output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-od", output_directory, PATH_MAX, &is_set));
    PetscCall(PetscLogStageRegister("Assembly", &stage_assembly));
    PetscCall(PetscLogStageRegister("Timestepping", &stage_timestep));

    dt = end_time / timesteps;

    PetscCall(PetscPrin2i(PETSC_COMM_WORLD, "timesteps", &timesteps, 1));
    PetscCall(PetscPrin2f(PETSC_COMM_WORLD, "dt", &dt, 1));
    PetscCall(PetscPrin2f(PETSC_COMM_WORLD, "end time", &end_time, 1));

    HDG hdg((1 << iteration) * space_dim, {tau, theta, dt});
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "true");
    hdg.plot_option("scale", "0.95");

    temp = hdg.make_initial(hdg.zero_vector());
    N = temp.size();
    hdg.plot_solution(temp, 0.); // needs petsc

    PetscCall(VecCreate(PETSC_COMM_SELF, &sol));
    PetscCall(VecSetSizes(sol, PETSC_DECIDE, N));
    PetscCall(VecSetType(sol, VECSEQ));

    PetscLogStagePush(stage_assembly);
    sparse_mat<std::vector<PetscReal>> mat_raw = hdg.trace_to_flux_mat(0.); // needs petsc
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 1, NULL, &mat));
    PetscCall(MatSetPreallocationCOO(mat, mat_raw.value_vec.size(), (PetscInt*)mat_raw.row_vec.data(), (PetscInt*)mat_raw.col_vec.data()));
    PetscCall(MatSetValuesCOO(mat, mat_raw.value_vec.data(), INSERT_VALUES));
    PetscLogStagePop();

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(PetscPrintf(PETSC_COMM_SELF, "timestepping...\n"));
    PetscLogStagePush(stage_timestep);
    for (PetscInt i = 0; i < timesteps; i++) {
        temp = hdg.residual_flux(hdg.zero_vector(), (i+1)*dt); // needs petsc
        PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, temp.size(), temp.data(), &rhs));
        PetscCall(VecScale(rhs, -1.));
        PetscCall(KSPSolve(ksp, rhs, sol));
        PetscCall(VecDestroy(&rhs));

        PetscCall(KSPGetIterationNumber(ksp, &its));
        iterations += its;

        PetscReal* sol_arr;
        PetscCall(VecGetArray(sol, &sol_arr));
        std::copy(sol_arr, sol_arr+N, temp.begin());
        PetscCall(VecRestoreArray(sol, &sol_arr));

        hdg.set_data(temp, (i+1)*dt);
        hdg.plot_solution(temp, (i+1)*dt); // needs petsc
    }
    PetscLogStagePop();

    temp2 = hdg.errors(temp, end_time);
    temp3 = hdg.norms(temp, end_time);
    for (size_t i = 0; i < temp3.size(); i++)
      temp3[i] = temp2[i] / temp3[i];
    avg_it = ((PetscReal)iterations) / timesteps;

    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final abs error", temp2.data(), temp2.size()));
    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "final rel error", temp3.data(), temp3.size()));
    PetscCall(PetscPrin2i(PETSC_COMM_SELF, "tot iterations", &iterations, 1)); 
    PetscCall(PetscPrin2f(PETSC_COMM_SELF, "avg iterations", &avg_it, 1)); 
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "wrote output to '%s/%s.*.vtu'\n", output_directory, output_filename));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}

    // char network_path[PATH_MAX];
    // PetscOptionsGetString(NULL, ns_pre, "network", network_path, PATH_MAX, &is_set);
