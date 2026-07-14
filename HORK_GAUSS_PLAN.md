# hoRK / Gauß time stepping for timowave — working plan

Goal: high-order (order-2s) Gauß time stepping for the condensed HDG Timoshenko
wave solver, per `hdg_gauss.pdf` (= hoRK `hdg.pdf`), **condense-last** so spatial
order stays `h^{p+1}` (condense-first / separate `Mbar,Kbar` caps at O(h²), cf.
the branch `jprecond/petscts-rewrite` `SecondOrder` loop).

## Settled decisions
- Build **fresh on ne9** (clean condense-last CN already lives here).
- Stepper is **hoRK** (`~/phd/hoRK/horkirk.c`), a custom PETSc `TS` type with
  `TSHORKIRKSetStageOperator(ts, fn, ctx)`, `fn(ts, sigma, &S, ctx)` returning the
  shifted stage operator for the **complex** `sigma = Δt·λ_i`. It hands us `sigma`
  and factors whatever `S` we return (KSPPREONLY+PCLU).
- horkirk's global state **u = the trace λ** (displacement/rotation, N dofs, as
  today). The dynamic `(y,z)=(u/r, v/s)` history lives in the per-edge local
  `data` (as today's CN driver already does), reconstructed each step.
- horkirk needs **complex PETSc** even for s=1 (its tableau builder calls `zgeev`;
  author's own note says so). ne9 is real today.

## Sequencing (real-first de-risk)
Two independent risks: (1) the novel HDG condense-last mapping, (2) the complex
build + horkirk plumbing. Do (1) first, in the existing **real** build, with a
hand-rolled s=1 step (no horkirk), then (2).

- **Phase A (real, no horkirk):** implement the per-stage condensed operator
  `Â(h)` and the field-history stage RHS `F̂`; drive **s=1 implicit midpoint**
  (`h = Δt/2`, `σ₁ = 2/Δt`) by hand; confirm 2nd-order temporal + `h^{p+1}` spatial
  convergence, i.e. the same behaviour as today's CN (ne9-01/02).
- **Phase B (complex + horkirk):** complex PETSc build; instantiate the local
  solver on `PetscComplex` (LAPACK → `zgetrf`); copy `horkirk.c` in; wire the
  stage-operator callback = `Â(σ)`, the RHS via IFunction, a TS post-step to
  reconstruct `(y,z)`. Verify order 2s for s=2,3.

## Key derivation — the stage operator `Â(h)`
The ne9 CN condensed displacement-Schur operator (`assemble_schur`, `S_entry`) is

    S = θτF + θ·C_sig·K + (C_u/(θΔt²)) M
      = θ · [ τF + C_sig·K + (C_u/h²) M ],     h := θΔt,   K=(B−G)M⁻¹G

so, dropping the irrelevant overall real scale θ,

    Â(h) = τF + C_sig·K + (C_u / h²) M.

This matches the note's per-stage form eq (12), `c(Y,Y) + (Δtθ_ℓ)²[a(Q,Q)+τ⟨…⟩]`,
after factoring `(Δtθ_ℓ)²`, with **h = Δt·θ_ℓ** (θ_ℓ = Butcher eigenvalue = horkirk's
λ_i; horkirk's `sigma = Δt·λ_i = h`). At CN, θ_ne9 = ½ ⇒ h = Δt/2 = Δt·θ_1. ✓

Coupled (shear/bend) blocks scale identically (`fill_coupled`: the θ factors out
of every entry), so the whole condensed operator is `θ·Â(h)`. **Â(h) is the single
knob**: real `Δt/2` at CN, complex `Δt·θ_ℓ` for s≥2. Note it is *not* affine in h
(the `1/h²` mass term) — condensation of the field-level affine `M_field−h·C_field`
makes the trace operator rational; that is expected and correct (Prop 1).

## Key reformulation — the stage RHS `F̂`
Today's `residual_flux`/`assemble_rhs_from_global_rhs` build the **θ-method**
(trapezoidal) load: θ/(1−θ) weighting of loads + accumulated `flux_*` + `v_old/Δt`.
The stage form (note eqs 7,11) instead takes history only through the initial
`(yⁿ,zⁿ)`:

    ĝ_y(w_y) = (f̃, w_y) + σ_ℓ ω_ℓ (zⁿ, w_y)
    ĝ_z(w_z) =            − σ_ℓ ω_ℓ (yⁿ, w_z),     σ_ℓ = 1/(Δt θ_ℓ) = 1/h,  ω_ℓ = (T⁻¹𝟙)_ℓ

then `F̂` is the condensed rhs of the local stage solve (5) with ζ=0. For s=1
(ω₁=1, σ₁=2/Δt) this is the implicit-midpoint load — 2nd order, differs from the
current trapezoidal load only at O(Δt²) in the forcing (both converge to the CN
result). Reconstruct `(yⁿ⁺¹,zⁿ⁺¹)` from the stage solution and store per edge
(reformulated `set_data`).

## Why not stock PETSc `TSTHETA` for Phase A
`TSTHETA` assumes an affine Jacobian `shift·M − C` with fixed `M,C` and reuses that
`M` in its `M·u̇ − C·u − g` residual. The condense-last `Â(h)` is **non-affine** in
the shift — the local solvers depend on `h` (they invert mass *and* stiffness
together, which is what buys `h^{p+1}`), so no fixed `(M,C)` gives `Â(h)=shift·M−C`
(it's `~1/h²`, and stays non-affine even in `α=h²`). Stock integrators can only do
the condense-**first** version (fixed `Mbar,Kbar` → the branch driver → O(h²)).
Re-assembling the operator per (complex) shift is exactly what horkirk's
`SetStageOperator` is for; `TSTHETA` has no equivalent. So Phase A is hand-rolled
(or reuses today's ne9 θ-step as the condense-last reference).

## Phase A code steps (real build, verifiable each step)
1. `Â(h)`: generalize `assemble_schur` to take the stage factor `h` instead of
   `(θ_,Δt_)`; assert it reproduces today's operator at `h=Δt/2` (compare against
   `trace_to_flux_mat`, or the `-loc_lu_full` full-matrix path).
2. `F̂`: stage RHS from `(yⁿ,zⁿ)` per eqs (7),(11); s=1 midpoint.
3. reconstruct/store `(y,z)` from the stage solution.
4. hand-rolled s=1 driver mode; run a conv-t/conv-x-style test; confirm 2nd-order
   temporal + `h^{p+1}` spatial (compare shape to `output/GOLDEN-ne9-0*`).

Verification helper: `experiments/ne9-conv-metrics.sh <results.json>`.
