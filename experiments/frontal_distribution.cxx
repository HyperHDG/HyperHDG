#include <petsc.h>
#include <cholmod.h>

#include <cmath>
#include <cstdio>
#include <vector>

static const char help_msg[] =
    "compute the distribution of supernodal frontal sizes, weighted by flops\n"
    "  -f <file>   : PETSc binary matrix (default mat.bin)\n"
    "  -csv <file> : write per-supernode data (cols,rows,flops) as csv\n";

static const char *ordering_name(int ordering) {
  switch (ordering) {
    case CHOLMOD_NATURAL: return "natural";
    case CHOLMOD_GIVEN: return "given";
    case CHOLMOD_AMD: return "amd";
    case CHOLMOD_METIS: return "metis";
    case CHOLMOD_NESDIS: return "nesdis";
    case CHOLMOD_COLAMD: return "colamd";
    case CHOLMOD_POSTORDERED: return "postordered";
    default: return "?";
  }
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, help_msg));

  char file[PETSC_MAX_PATH_LEN] = "mat.bin";
  char csv_file[PETSC_MAX_PATH_LEN] = "";
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", file, sizeof(file), NULL));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-csv", csv_file, sizeof(csv_file), NULL));

  Mat A;
  PetscViewer viewer;
  PetscCall(MatCreate(PETSC_COMM_SELF, &A));
  PetscCall(MatSetType(A, MATSEQAIJ));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, file, FILE_MODE_READ, &viewer));
  PetscCall(MatLoad(A, viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  const PetscInt *ia, *ja;
  PetscInt n;
  PetscBool done;
  PetscCall(MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done));
  PetscCheck(done, PETSC_COMM_SELF, PETSC_ERR_SUP, "MatGetRowIJ gave no ij structure");
  static_assert(sizeof(PetscInt) == sizeof(int32_t), "CHOLMOD_INT needs 32-bit PetscInt");

  cholmod_common c;
  cholmod_start(&c);
  c.supernodal = CHOLMOD_SUPERNODAL;  // force supernodal so L carries super/pi

  // CSR of the (structurally symmetric) matrix reused as CSC; stype=1 keeps the
  // upper triangle and ignores the rest, so no explicit transpose is needed.
  cholmod_sparse A_ch = {};
  A_ch.nrow = A_ch.ncol = (size_t)n;
  A_ch.nzmax = (size_t)ia[n];
  A_ch.p = (void *)ia;
  A_ch.i = (void *)ja;
  A_ch.stype = 1;
  A_ch.itype = CHOLMOD_INT;
  A_ch.xtype = CHOLMOD_PATTERN;
  A_ch.dtype = CHOLMOD_DOUBLE;
  A_ch.sorted = 1;
  A_ch.packed = 1;

  cholmod_factor *L = cholmod_analyze(&A_ch, &c);
  PetscCheck(L && c.status == CHOLMOD_OK, PETSC_COMM_SELF, PETSC_ERR_LIB,
             "cholmod_analyze failed (status %d)", c.status);
  PetscCheck(L->is_super, PETSC_COMM_SELF, PETSC_ERR_LIB, "analysis is not supernodal");

  const int32_t *super = (const int32_t *)L->super;
  const int32_t *pi = (const int32_t *)L->pi;
  const size_t nsuper = L->nsuper;

  FILE *csv = NULL;
  if (csv_file[0]) {
    csv = fopen(csv_file, "w");
    PetscCheck(csv, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "cannot open %s", csv_file);
    fprintf(csv, "cols,rows,flops\n");
  }

  // log2 bins over the front height (rows = leading dimension of the dense panel)
  constexpr int NBINS = 32;
  double bin_flops[NBINS] = {};
  long bin_count[NBINS] = {};
  double total_flops = 0;
  long max_rows = 0, max_cols = 0;

  for (size_t s = 0; s < nsuper; ++s) {
    const long cols = super[s + 1] - super[s];
    const long rows = pi[s + 1] - pi[s];
    // multiply-add pairs of the dense panel factorization: column k of the
    // supernode touches rows-k rows, so potrf+trsm+syrk sum to sum_k (rows-k)^2
    double flops = 0;
    for (long k = 0; k < cols; ++k) flops += (double)(rows - k) * (rows - k);
    if (csv) fprintf(csv, "%ld,%ld,%.0f\n", cols, rows, flops);

    int b = 0;
    while ((1L << (b + 1)) <= rows && b < NBINS - 1) ++b;
    bin_flops[b] += flops;
    bin_count[b] += 1;
    total_flops += flops;
    if (rows > max_rows) max_rows = rows;
    if (cols > max_cols) max_cols = cols;
  }
  if (csv) fclose(csv);

  PetscCall(PetscPrintf(PETSC_COMM_SELF,
                        "n=%" PetscInt_FMT "  nnz=%" PetscInt_FMT "  ordering=%s  nsuper=%zu\n"
                        "lnz=%.3e  flops(panels)=%.3e  flops(cholmod fl)=%.3e\n"
                        "max front: rows=%ld cols=%ld\n\n",
                        n, ia[n], ordering_name(L->ordering), nsuper, c.lnz, total_flops, c.fl,
                        max_rows, max_cols));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "%18s %12s %14s %8s %8s\n", "front rows", "supernodes",
                        "flops", "%flops", "cum%"));
  double cum = 0;
  for (int b = 0; b < NBINS; ++b) {
    if (!bin_count[b]) continue;
    cum += bin_flops[b];
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "[%7ld, %7ld) %12ld %14.3e %8.2f %8.2f\n", 1L << b,
                          1L << (b + 1), bin_count[b], bin_flops[b], 100 * bin_flops[b] / total_flops,
                          100 * cum / total_flops));
  }

  cholmod_free_factor(&L, &c);
  cholmod_finish(&c);
  PetscCall(MatRestoreRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
