#include <stdio.h>
#include <petsc.h>

#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include "parameters.hxx"

static const char help[] = "experiments regarding timoshenko networks\n";

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
    int errcode = 0;

    constexpr unsigned int poly_deg = 5;
    using Top = Topology::File<1,3>;
    using Geo = Geometry::File<1,3>;
    using NDes = NodeDescriptor::File<1,3>;
    // using LSol = LocalSolver::TimoshenkoBeam<1,3,poly_deg,2*poly_deg,LocalSolver::TimoschenkoBeamParametersClamped>;
    using LSol = LocalSolver::Diffusion<1,5,10,ConstantDiffusionParameters>;
    using HDG = GlobalLoop::Elliptic<Top,Geo,NDes,LSol>;
    constexpr unsigned int n_dofs_per_node = 6;

    char output_directory[PATH_MAX] = "output";
    char output_filename[PATH_MAX] = "network";
    char domain_filepath[PATH_MAX] = "domains/grid_2.geo.bin.zstd";

    PetscLogStage s_as, s_it, s_rf;

    PetscBool is_set;
    PetscInt N;
    PetscReal err, sol_norm;
    PetscInt iterations;

    std::vector<PetscReal> temp, temp2, temp3, zero_v;
    sparse_mat<std::vector<PetscReal>> mat_coo;
    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "initialization...\n"));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-o", output_filename, PATH_MAX, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-od", output_directory, PATH_MAX, &is_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-domain", domain_filepath, PATH_MAX, &is_set));
    PetscCall(PetscLogStageRegister("Assembly", &s_as));
    PetscCall(PetscLogStageRegister("Iteration", &s_it));
    PetscCall(PetscLogStageRegister("residual_flux", &s_rf));

    HDG hdg(domain_filepath);
    hdg.plot_option("fileName", output_filename);
    hdg.plot_option("outputDir", output_directory);
    hdg.plot_option("printFileNumber", "false");
    hdg.plot_option("scale", "1");

    zero_v = hdg.zero_vector();
    N = zero_v.size();
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &sol));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &rhs));

    PetscLogStagePush(s_as);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "assembly...\n"));
    mat_coo = hdg.trace_to_flux_mat();
    PetscCall(MatCreateSeqAIJFromTriple(PETSC_COMM_SELF, N, N,
      (PetscInt*)mat_coo.row_vec.data(), (PetscInt*)mat_coo.col_vec.data(),
      mat_coo.value_vec.data(), &mat, mat_coo.value_vec.size(), PETSC_FALSE /* 0-based */));
    PetscLogStagePop();

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, mat, mat));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE)); // no diagonal preconditioning
    PetscCall(KSPSetFromOptions(ksp));

    PetscCall(PetscPrintf(PETSC_COMM_SELF, "iteration...\n"));
    std::span<PetscReal> rhs_span;
    std::span<PetscReal> sol_span;
    PetscCall(VecGetSpan(rhs, rhs_span));
    PetscCall(VecGetSpan(sol, sol_span));

    PetscLogStagePush(s_rf);
    hdg.residual_flux2(zero_v, rhs_span, 0.);
    PetscLogStagePop();
    PetscCall(VecScale(rhs, -1.));

    PetscLogStagePush(s_it);
    PetscCall(KSPSolve(ksp, rhs, sol));
    PetscLogStagePop();

    PetscCall(KSPGetIterationNumber(ksp, &iterations));

    hdg.plot_solution(sol_span);

    PetscCall(VecRestoreSpan(rhs, rhs_span));
    PetscCall(VecRestoreSpan(sol, sol_span));

    PetscCall(PetscPrin2i(PETSC_COMM_SELF, "iterations", &iterations, 1)); 
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "wrote output to '%s/%s.*.vtu'\n", output_directory, output_filename));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&mat));
    PetscCall(VecDestroy(&sol));

    PetscCall(PetscFinalize());
    return 0;
}

    // char network_path[PATH_MAX];
    // PetscOptionsGetString(NULL, ns_pre, "network", network_path, PATH_MAX, &is_set);
