# Gauß time stepping for timowave — plan v2 (no horkirk, no PETSc TS)

Goal: order-2s fully implicit Gauß collocation for the condensed HDG Timoshenko
wave solver, condense-last so spatial order stays `h^{p+1}`, per `hdg_gauss.pdf`
(= `hdg.pdf`). v1's horkirk / PETSc-TS route is dropped — conforming to a TS
interface bought only plumbing; the stage loop is hand-rolled in the driver.

## Why the loads simplify (settled understanding)

Gauß is *collocation*: the saddle system is imposed exactly at interior nodes,
never at `t^n`. The spatial operator is only ever applied to stage unknowns;
history enters solely through the time derivative of the collocation
polynomial, i.e. through mass inner products of `(y^n, z^n)` (hdg.pdf eq 7):

    ĝ_y(w) = (f̃_ℓ, w) + σ_ℓ ω_ℓ (z^n, w),   ĝ_z(w) = −σ_ℓ ω_ℓ (y^n, w),
    σ_ℓ = 1/(Δt θ_ℓ)

Contrast CN/trapezoid (main.pdf p.15): the old state appears *under the
operator* (`b(q^{k−1},·)`, `τ⟨y^{k−1}−λ^{k−1},·⟩`, `d(z^{k−1},·)`), which is
what dragged `q_old`/`λ_old` (the `flux_*` caches) into `data`. `q=(n,m)` and
the trace `λ` are algebraic — slaved to the instantaneous state, recomputed per
stage, never history. The `flux_*` bookkeeping is trapezoid-specific and does
not carry over.

## Architecture (settled)

- **One global loop, one geometry, one data container.** No per-stage loop
  instances, no state sharing between siblings.
- **Stage index as a defaulted extra parameter** (`stage = 0`) on the local
  solver's `trace_to_flux` / `residual_flux` / `set_data` and the loop's
  `trace_to_flux_mat` / residual / set_data entries. No `set_active_stage`
  member state: solver stays `const`, stage matrices assemble in parallel.
- **`n_stages` compile-time** template parameter (like `poly_deg`). Tableau
  data are *solver members* (not `data_type`), computed once in the
  constructor: hardcoded Gauß Butcher `A, b` for s = 1..3, eigen-decomposition
  `A = TΘT⁻¹` via LAPACK `dgeev`, conjugate eigenvector columns paired.
- **`theta` keeps name and meaning**: `θ_ℓ` = Butcher eigenvalue; stage shift
  `h_ℓ = Δt·θ_ℓ`; s=1 gives `θ₁ = 1/2` — today's `theta_ = .5` semantics.
- **Only ⌈s/2⌉ representatives** (convention `Im θ_ℓ ≥ 0`) are solved and
  stored. Conjugate partners are analytic (`ζ̄, ȳ, z̄` since A, forms, data are
  real) and never materialized.
- **Per-representative tableau quantities**: `θ_ℓ`; `ω_ℓ = (T⁻¹𝟙)_ℓ` (loads);
  recombination weight `w_ℓ = (dᵀT)_ℓ` with `dᵀ = bᵀA⁻¹`; multiplicity
  `m_ℓ ∈ {1,2}`. Endpoint update:
      y⁺ = (1 − Σ_j d_j)·y^n + Σ_ℓ m_ℓ·Re(w_ℓ · y_ℓ)       (same for z)
  s=1 sanity: d = 2, w₁ = 2 → y⁺ = 2y₁ − y^n (midpoint extrapolation).
- **`data_type`**: real state `(u,r,v,s)` + per-representative complex LU
  caches + per-representative complex stage-field slots `(y_ℓ, z_ℓ)`.
  POD, fixed size via `n_stages` (redistribution-safe).
- **Stage operator — LANDED, single path** (commits 404feab4, 1a7a121e): the
  whole local solver is assembled in stage normalization by ONE
  `assemble_loc_matrix<with_z>(sigma)` (`with_z=false` = the initializers'
  static system); no flags, no `h` — the step enters only through
  `sigma = 1/(theta*Δt)`. Deviation from hdg.pdf kept deliberately: z-rows
  carry the C_u scaling ((z,z)=M, (z,y)=−σ·C_u·M) so massless welds (C_u=0)
  stay regular (z≡0); costs the literal saddle symmetry, not the condensed
  operator. Schur path (`assemble_schur`: S = τF + C_sig·K + σ²C_u·M) and
  full-LU path both verified against GOLDEN (θ=.5, θ=1, -loc_lu_full) to
  printed precision. The trapezoid loads are row-rescaled into this
  normalization at the end of `assemble_rhs_from_global_rhs` — that block is
  what the Gauss loads (eq 7) replace. `recover_dual(y, λ)` (algebraic n,m
  recovery) and `add_dirichlet_rhs_static` exist and are shared by the unified
  initializers; recover_dual is the output-time diagnostics helper too.
- **Step protocol** (driver-controlled, no tally):
  1. per representative ℓ: `b_ℓ = residual_flux(stage ℓ)`;
     solve `Â_ℓ ζ_ℓ = b_ℓ` — independent, any order, parallelizable;
  2. `set_data(ζ_ℓ, ℓ)`: solve stage locals from ζ_ℓ + stored `(y^n,z^n)`,
     stash `(y_ℓ, z_ℓ)` in slot ℓ. Idempotent (re-call overwrites).
  3. `finalize_step()`: pure real axpys over slots → new `(y,z)`. Explicit
     call — completion is decided by the driver, not inferred per edge.
  Debug `hy_assert`: imaginary residue of the finalize combination ≈ 0
  (catches any conjugation-convention slip in the θ/ω/w chain).

## Landed: Gauss s=1 REPLACES the trapezoid scheme (no dual path)

The trapezoidal θ-scheme is ripped out, not kept alongside: no `-gauss` /
`-theta` options, no `flux_*` caches, no `compute_fluxes`, no θ-averaged
load integrators, no `assemble_rhs_from_global_rhs`. The single time stepper
is the Gauss step protocol above; `timowave.cxx` runs it unconditionally
(stage solve → `set_data` stash → `finalize_step` → endpoint trace
`λ⁺ = 2ζ₁ − λⁿ`, exact at s=1 only).

- `n_stages` is a compile-time template parameter of `TimoshenkoWave`,
  tied to the spatial degree in the driver alias: deg ≤ 2 → s = 1, deg 3 →
  s = 2 (temporal order 2s covers spatial order p+1).
- Tableau data are built at construction by `Gauss::build_tableau`
  (`gauss_tableau.hxx`, ported from `~/phd/hoRK/horkirk.c`): Gauss nodes by
  Newton on P_s, Butcher A = W V⁻¹, zgeev eigen-decomposition with enforced
  conjugate pairing and phase-fixed real columns. Value-form weights:
  ω = T⁻¹𝟙, w = (bᵀA⁻¹T), affine = 1 − bᵀA⁻¹𝟙 = R(∞) = (−1)^s.
  Invariants tested in `tests_c++/gauss_tableau.cxx`.
- Verified s=1 (ne9-06/ne9-07 scripts, output saved): temporal order 2 and
  spatial h^{p+1} against the θ=0.5 rows of `output/GOLDEN-ne9-0*`;
  iteration counts identical (same operator), e_rel differs at O(Δt²) only
  (midpoint vs endpoint-averaged sampling of time-dependent Dirichlet data).

## Landed: s = 2 (order 4) via complex stage solves

- Stage scalar `stage_float_t` = complex for s ≥ 2; only the ⌈s/2⌉
  representatives are solved (conjugates analytic). Local stage systems:
  per-representative complex full-matrix LU cached in `data_type`
  (`[[no_unique_address]]`-gated so s = 1 keeps its Schur/full-LU paths and
  data size). Stage loads combine the data at ALL collocation nodes with the
  representative's T⁻¹ row (eq 7): rhs_ℓ = Σⱼ τ_{ℓj}·data(tⁿ+cⱼΔt) +
  σ_ℓω_ℓ·(history mass terms).
- NO stage-specific global-loop plumbing: the stage index rides inside a
  `Gauss::StageTime{time, stage}` passed as the (now template-typed) `time`
  argument of the GENERIC loop entries (`trace_to_flux_mat` /
  `residual_flux2` / `set_data`), instantiated with complex vectors; the
  local solver unpacks it via `split_stage_time`. The prototype macros
  derive the scalar from the span/matrix type.
- Global stage solves: NO complex PETSc — the condensed complex stage
  operator is assembled once (complex COO through the generic probing),
  densified and zgetrf-factored; each step is one zgetrs (conv-study scale;
  production choice complex-PETSc vs 2N×2N real block deferred). PETSc real
  Vecs shuttle Re/Im halves for scatters/layout only. Single rank enforced.
- Endpoint trace: λ⁺ = affine·λⁿ + Σ mult·Re(w·ζ) (driver AXPYs via
  `stage_weights`). Same recombination as the state — Gauss-quadrature
  superconvergent; observed e_trace order ≈ 3.3 on wave4 (time-dependent
  Dirichlet), the static-trace-solve upgrade remains the open item.
- **Verified s=2** (`ne9-08-conv-gauss2.sh`, deg 3): temporal e_rel ratios
  ≈ 2^4.0 over nt 8→32, flooring exactly at the golden deg-3 spatial floor;
  spatial h⁴ values match `GOLDEN-ne9-02` deg-3 rows to 3–4 digits with
  nt = 64 instead of the golden's nt = 8000.

## Landed: complex PETSc build + one unified driver path

- Second spack env `spack/complex` (petsc+complex+mumps+…) with CMake preset
  `complex` (CMakeUserPresets.json, local like `openblas`; cache var
  `HYPERHDG_COMPLEX=ON` builds only prin2 + timowave):
      eval $(spack env activate --sh spack/complex) && cmake --preset complex
      cmake --build --preset complex --target timowave
- The driver has ONE time-stepping path for any s: per-representative stage
  operators (PETSc Mats from the complex COO of the generic probing) and one
  KSP each. s = 1 is real-valued and runs in both builds with the classic
  defaults (CG + net2as); s ≥ 2 requires the complex build and defaults to
  direct LU (the stage operators are complex symmetric, NOT Hermitian — CG /
  cholmod theory does not apply; net2as stays real-build-only for that reason
  plus ~10 mechanical PetscReal*-vs-PetscScalar* sites).
- The endpoint trace λ is protocol-managed per-edge state: set_data stashes
  the stage trace ζ_ℓ next to the stage locals, finalize_step recombines
  state AND trace; errors/energy read `data.lambda_old`, their span argument
  is vestigial. No real-valued trace crosses the driver boundary; the driver
  holds no λ vector. make_initial seeds λ⁰ (computed in a real buffer,
  widened into possibly-complex outputs). -static was removed from the
  driver/interface pending its rework.
- Complex↔real bridges where the s=1 real operator meets complex storage:
  `apply_local_flux_widened` (narrow → real machinery → widen) and the
  set_data narrowing branch; read_domain/plot/prin2 made scalar-safe (plot
  values are `plot_value_t = real`; legacy vtu writer gated for complex).
- All dense-LAPACK global-solve machinery, the Re/Im half-vectors, and the
  wrapper zipping glue are gone; `-mat_cache` removed.
- Verified: s=1 anchors bit-exact in BOTH builds (deg 2: e_rel/e_trace
  2.48402e-3 / 8.45959e-5); s=2 via PETSc complex LU identical to the former
  dense path (deg 3, nt 32: 1.81203e-5 / 4.70214e-5); ne9-08 script now runs
  the complex build.

## Next steps

1. Endpoint trace and (n,m) at output times via a static condensed trace
   solve + `recover_dual` (restores full-order e_trace for s ≥ 2 with
   time-dependent data).
2. s = 3 (order 6): tableau builder already generic; needs the real-eigenvalue
   representative wired through (mult = 1 slot) and a driver test.
3. Production-scale stage solves at network size: iterative solvers /
   preconditioning for the complex-symmetric stage operators (net2as-style DD
   is a research question there), np > 1 for the complex build.

Verification helpers: `experiments/ne9-conv-metrics.sh <results.json>`;
tableau invariants: `tests_c++/gauss_tableau.cxx`.
