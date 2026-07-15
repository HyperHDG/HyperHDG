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
- **Stage operator**: `Â(h) = τF + C_sig·K + (C_u/h²)M`, i.e. the θ-method
  operator with the overall `θ·(...)` normalization dropped consistently from
  matrix and rhs (ζ invariant to the common scale). At s=1 equals today's
  operator up to the overall factor θ = 1/2. **Landed** (merge of
  jprecond/hork-gauss): `assemble_loc_matrix_stage(h)` +
  `assemble_rhs_from_lambda_stage` behind `-loc_stage`, full-matrix LU form —
  chosen over a hand-rolled stage Schur because it complexifies trivially
  (zgetrf). Verified: symmetric to ~1e-17, CN run bit-identical to GOLDEN.
  A stage-Schur fast path is an optional later optimization.
- **Step protocol** (driver-controlled, no tally):
  1. per representative ℓ: `b_ℓ = residual_flux(stage ℓ)`;
     solve `Â_ℓ ζ_ℓ = b_ℓ` — independent, any order, parallelizable;
  2. `set_data(ζ_ℓ, ℓ)`: solve stage locals from ζ_ℓ + stored `(y^n,z^n)`,
     stash `(y_ℓ, z_ℓ)` in slot ℓ. Idempotent (re-call overwrites).
  3. `finalize_step()`: pure real axpys over slots → new `(y,z)`. Explicit
     call — completion is decided by the driver, not inferred per edge.
  Debug `hy_assert`: imaginary residue of the finalize combination ≈ 0
  (catches any conjugation-convention slip in the θ/ω/w chain).

## Phases (each ends in a commit + verification)

### Phase 0 — concepts refactor (independent cleanup, lands first)
Replace the `HAS_MEMBER_FUNCTION` trait layer (pre-`requires` era) with inline
C++20 `requires`-expressions at the call sites:
- `global_loop/prototype.hxx`: the three `prototype_*` macros lose their
  `has_fun_name` parameter; dispatch via
  `if constexpr (requires { local_solver_.fun_name(args...); })`.
- Direct `if constexpr` sites: parabolic, hyperbolic, elliptic,
  nonlinear_eigenvalue, shifted_inverse_eigenvalue, mass_approx_eigenvalue,
  plot.hxx (`bulk_values`/`energy`/`n_energy_components` — shared machinery,
  mechanical swap).
- Fallback branches become dependent `static_assert` (`always_false_v`)
  instead of Release-silent `hy_assert(false, …)` → the `error_def`
  silent-zero bug class becomes a compile error.
- Delete the macro from `compile_time_tricks.hxx`.
Verify: full build + ctest; rerun ne9-01/02 → `ne9-conv-metrics.sh` identical
to `output/GOLDEN-ne9-0*` (pure refactor).

### Phase 1 — complex LAPACK plumbing
`lapack_factorize` / `lapack_solve_factored` overloads for
`std::complex<double>` (`zgetrf`/`zgetrs`); `dgeev` wrapper for tableau setup.

### Phase 2 — stage machinery in timowave (real, s=1)
- Stage operator Â(h): DONE (`assemble_loc_matrix_stage` behind `-loc_stage`,
  see above).
- `n_stages` template param (default 1); tableau members; slots in `data_type`.
- Stage entry points alongside the untouched θ-path: `trace_to_flux(…,stage)`,
  `residual_flux(…,stage)` with eq-(7) loads, `set_data(ζ,stage)`,
  `finalize_step`.
- s=1 is real ⇒ end-to-end testable with zero complex infrastructure.
- **Verify — CN reproduction against GOLDEN:**
  (i) stage matrix at θ=1/2 equals today's `trace_to_flux_mat` up to the
  overall θ scale;
  (ii) rerun `ne9-01-conv-t` / `ne9-02-conv-x` through the stage path,
  compare `experiments/ne9-conv-metrics.sh` vs `output/GOLDEN-ne9-0*`:
  identical orders and floors. For f = 0 with time-independent BCs the
  midpoint and trapezoid update maps coincide → must match to solver
  tolerance; configs with time-dependent Dirichlet data may differ at
  O(Δt²) in the load (midpoint samples `t^{n−1/2}`, trapezoid averages
  endpoints) — same order, same floor.

### Phase 3 — global loop + driver stage plumbing
- `global_loop/hyperbolic.hxx`: thread `stage` through matrix / residual /
  set_data entries, add `finalize_step` pass; `requires`-dispatch with
  `static_assert` fallback for the new signatures.
- `experiments/timowave.cxx`: stage loop, per-stage factorized solves.

### Phase 4 — complex end-to-end (s = 2, 3)
- Loop instance with `dof_value_t = std::complex<double>` for stage solves;
  state and diagnostics stay real.
- Global stage solves: factor-once complex direct solve; conv studies do not
  need complex PETSc. Production choice (complex PETSc build vs 2N×2N real
  block form — latter loses cholmod/SPD) deferred until network-scale runs.
- Verify: temporal order 4 (s=2) / 6 (s=3) on ne9-01-style runs; spatial
  `h^{p+1}` unchanged (ne9-02 style); energy conservation over long runs.

## Open items
- **λ and (n,m) at output times**: algebraic → recompute on demand. Endpoint
  trace via a static condensed trace solve given the new `(y,z)` (do NOT
  extrapolate the collocation polynomial — loses accuracy at low s), then
  local recovery of `(n,m)`. Needed for `e_trace` and energy diagnostics;
  volume errors need only `(y,z)`.
- **Retire the trapezoidal path** (`flux_*`, `(n,m)`-as-state, θ-weighted
  rhs) once the Gauß path is the default; consolidate `data` state to the
  single real `(y,z)` block.

Verification helper: `experiments/ne9-conv-metrics.sh <results.json>`.
