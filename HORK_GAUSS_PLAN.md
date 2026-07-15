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

- `n_stages` is a compile-time template parameter of `TimoshenkoWave`
  (default 1; `static_assert(n_stages == 1)` until the complex plumbing).
- Tableau data are per-stage `std::array`s (`stage_theta_`, `stage_c_`,
  `stage_omega_`, `stage_w_` + scalar `stage_affine_`), currently the single
  implicit-midpoint entry {1/2, 1/2, 1, 2; −1} whose trajectory is CN's.
  s ≥ 2 fills them from the Butcher eigen-decomposition (dgeev) instead.
- `sigma(stage)`, `assemble_rhs_stage(edge, time, stage)`,
  `set_data(λ, edge, time, stage)`, per-stage `data.stage_coeffs[stage]`.
- Verified (ne9-06/ne9-07 scripts, output saved): temporal order 2 and
  spatial h^{p+1} against the θ=0.5 rows of `output/GOLDEN-ne9-0*`;
  iteration counts identical (same operator), e_rel differs at O(Δt²) only
  (midpoint vs endpoint-averaged sampling of time-dependent Dirichlet data).

## Next steps

1. Complex LAPACK plumbing: `zgetrf`/`zgetrs` overloads for the local stage
   solves, `dgeev` wrapper for the tableau setup.
2. s = 2, 3 tableaux + complex per-representative stage solves (local caches
   and `stage_coeffs` become complex, ⌈s/2⌉ representatives).
3. Global loop / driver stage plumbing: stage-indexed `trace_to_flux_mat` /
   residual entries, per-stage factorized global solves.
4. Endpoint trace and (n,m) at output times via a static condensed trace
   solve + `recover_dual` (the s=1 extrapolation shortcut does not carry
   over).

Verification helper: `experiments/ne9-conv-metrics.sh <results.json>`.
