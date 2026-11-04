# Local operator / global loop redesign, v2 — interfaces only

Status: design note, 2026-10-28; the prototype on branch operator-loop implements it (see python/README.md).
Decisions folded in from review: kind/args/key split, constraints in the
loop, block apply instead of matrix(), reduce in the loop only, parallel
prepare, optional update/operator_key, implicit store policy
(prepare/release vs streaming), SmallMat<n_loc, k> instead of a block
type.

## 1. Goals

- One global loop. Elliptic, Parabolic, MassApprox, ShiftedInverse,
  Hyperbolic, Generic become one scatter–gather engine.
- The local solver exposes operations, never block matrices. A, B, C, D,
  Schur reductions, diagonal mass blocks, complex stage operators stay
  private to the solver.
- Caching is a policy of the loop. The solver says how to precompute and
  what the precomputation depends on; the loop says when and for whom.
- Constraints (Dirichlet) are a runtime index set owned by the loop. The
  solver never masks.
- Every operation the loop may run in parallel is `const` and free of
  hidden mutation. Mutation happens only in `prepare` and `update`, both
  scheduled by the loop with exclusive access. This rule is what makes
  shared, distributed, and hybrid parallelism possible without touching
  the local interface.
- Minimal solver: one `local_operator` overload returning an object with
  `apply` and `residual`. Everything else is optional.

## 2. Vocabulary

    Kind   compile-time tag selecting which operator is meant
           struct Stiffness {}; struct Mass {}; struct Shifted {}; struct Stage {};
    Args   small runtime struct passed to every call
           { double time; std::complex<double> sigma; int stage; const Trace* lin; ... }
           Each solver defines the Args it understands; the loop forwards it.
    Key    the subset of (edge class, Args) the precomputed data depends on.
           Returned by the solver; compared by the loop.
    State  mutable per-edge data advanced by the time stepper
           (old coefficients, old trace, stage stashes). Lives in
           hyper_edge.data. Never read by the loop.
    n_loc  2 * hyEdge_dim * n_dofs_per_node, the trace size of one edge.
    Trace  SmallVec<n_loc, value_t>; face/dof are i / n_dofs and i % n_dofs.
    SmallMat<n_loc, k, value_t>  k traces as columns (column-major, as in
           dense_la). k = 1 is a Trace; k = n_loc with the identity is
           assembly and its output is the dense local matrix. No new
           types: both already live in dense_la.hxx.

## 3. LocalOperator (immutable, owned by the loop's store)

A value returned by the solver. Holds what the solver wants to precompute
for repeated application: LU factors, a Schur reduction, tabulated
boundary integrals, or nothing. It is a pure function of (geometry, node
descriptor, solver constants, Key). It never contains State and never
depends on the parts of Args outside its Key, so one instance can be
shared across equivalent edges and across time steps.

Required members:

    // y = J δ for k columns at once. For linear problems J = S.
    // k = 1 is the matrix-free product; k = n_loc with the identity block
    // is assembly. The solver may hoist boundary integrals across columns
    // and use a multi-rhs back-substitution; a naive solver loops over
    // columns. One implementation serves both paths, so the assembled
    // matrix equals the matrix-free operator by construction.
    template <unsigned k>
    void apply(const Edge&, const SmallMat<n_loc, k, value_t>& delta,
               const Args&, SmallMat<n_loc, k, value_t>& y) const;

    // r = F(λ): residual of the condensed problem on this edge. Reads
    // State through `edge` (history) and Args (time-dependent loads).
    void residual(const Edge&, const Trace& lambda, const Args&, Trace& r) const;

Optional members (checked with `requires`):

    // Local solution for a given trace; needed only by reductions that
    // evaluate volume quantities (e.g. Diffusion L2 error). Solvers that
    // keep coefficients in State do not provide it.
    SmallVec<n_loc_dofs, value_t>
    local_solution(const Edge&, const Trace& lambda, const Args&) const;

No matrix(), no constrained(). Frame changes (node↔edge dofs), τ-terms,
Schur back-substitution happen inside apply/residual. Dirichlet masking
does not happen here at all (see §5.3).

Why apply + residual and not A, B, C, D:
- mirrors the loop, which only needs F and J·δ,
- covers linear, affine (time stepping), and nonlinear local problems,
- lets each solver pick its own decomposition and caching.

Why block apply and not matrix():
- assembly is apply with the identity block; one code path,
- the block form keeps the integral hoisting of the hand-written
  Diffusion::trace_to_flux_mat and adds a level-3 multi-rhs solve,
- block Krylov, several stage vectors, several eigenvectors reuse it.

## 4. LocalSolver (stateless apart from edge-independent constants)

Required:

    static constexpr unsigned n_dofs_per_node;
    using State = ...;                       // may be empty
    using Args  = ...;                       // the runtime parameter struct

    // Build the operator for one edge and one Kind. Only place where
    // factorizations happen. Called from prepare, never from apply/residual.
    // One overload per supported Kind; a missing Kind is a compile error
    // at the loop's call site.
    LocalOperator local_operator(const Edge&, Stiffness, const Args&) const;

Optional (checked with `requires`, defaults in the loop):

    // Which parts of (edge, Args) the operator depends on. Default: the
    // edge index, i.e. one operator per edge, no sharing. Examples:
    //   Diffusion, constant coefficients, UnitCube : geometry class only
    //   Shifted                                    : (geometry class, sigma)
    //   TimoshenkoWave                             : (edge, stage), not time
    Key operator_key(const Edge&, Kind, const Args&) const;

    // Edge-exclusive state transition (set_data / finalize_step today).
    // Default: absent; GlobalLoop::update is then a compile error if called.
    void update(Edge&, const Trace& lambda, const Args&) const;

Per-edge quantities for reductions (errors, energies, plot values) are
ordinary named members with whatever signature suits them; they are
called from functors passed to GlobalLoop::reduce (§5.5), so they are not
part of the interface.

Minimal solver example (steady, linear, one Kind, no sharing, no state):

    struct Args { double time = 0.; };
    struct State {};
    LocalOperator local_operator(const Edge& e, Stiffness, const Args& a) const;
    // LocalOperator provides apply<k>() and residual().

## 5. GlobalLoop

    template <class HyperGraph, class LocalSolver>
    class GlobalLoop {
      // y = J x, matrix-free.
      template <class Kind> void apply(Kind, const Args&, const Vec& x, Vec& y);

      // r = F(x).
      template <class Kind> void residual(Kind, const Args&, const Vec& x, Vec& r);

      // Triplets of J, via apply<n_loc> with the identity block.
      template <class Kind> Triplets assemble(Kind, const Args&);

      // State transition for all edges; the only phase mutating State.
      // Requires LocalSolver::update.
      void update(const Args&, const Vec& x);

      // Global reductions, see §5.5.
      template <class Kind, class F, class Combine = plus>
      auto reduce(Kind, const Args&, const Vec& x, F&& f, Combine&& c = {});

      // Store control, see §5.2. Optional: without prepare, sweeps stream.
      // After prepare, the caller owns the store and must call release.
      template <class Kind> void prepare(Kind, const Args&);
      void release();

      // Constraints, see §5.3.
      void set_constraints(ConstraintSet);
    };

### 5.1 Scatter–gather core

For every edge: gather its trace from x, look up its operator, call
apply/residual, scatter into y. The index mapping between edge dofs and
global dofs is owned by the loop and visible to it (needed for coloring,
distribution, assembly). Only the frame rotation inside an edge belongs
to the solver.

### 5.2 Store and prepare (implicit policy)

There is no policy class. The loop is either prepared for a Kind or not.

Prepared: the user called prepare(Kind, Args). Operators stay in the
store until release(), or until the next prepare of the same Kind, which
rebuilds that Kind's store from scratch (so a sweep over shifts does not
accumulate operators). Users who call prepare must call release when
done; documented at the function.

Not prepared: apply/residual/reduce/assemble stream. Each thread builds
the operator for its current edge, uses it, drops it. Memory O(1), no
map, no sharing by key. This is the default for one-off applications and
for simple users who never learn about prepare.

prepare(Kind, Args) in two passes:
  1. compute operator_key for every edge; deduplicate into (Key →
     representative edge, slot); store each edge's slot index.
  2. parallel over slots: build the operator from the representative
     edge. Pure, no shared writes.
The map is never touched inside a parallel region; apply reads slot
indices only. Under MPI all of this is rank-local.

Footgun: a Krylov loop calling apply repeatedly without prepare rebuilds
every time. It is slow, not wrong. Mitigations: the Python drivers that
run a solve call prepare/release themselves; in debug mode the loop warns
after N unprepared sweeps of the same Kind.

With the default key (edge index) a prepared store is per-edge caching,
which is timowave's current behaviour with the lazy fill made explicit.
With a geometry-derived key it is #105's factorize_once, correct by
construction, no sampling check.

### 5.3 Constraints

The loop owns a ConstraintSet: a vector of global dof indices, built once
at runtime from a user predicate over (node descriptor, dof component).
The ternary face digits in the geo file stay as data; which faces and
components are constrained is decided at runtime, e.g. "only the right
boundary in dimension 1".

Prescribed values live in the trace vector: the driver writes Dirichlet
data into the constrained entries before solving (and rewrites them each
step if time dependent). The local operator receives the full trace,
which is what the local problem wants.

After the scatter:
  apply     : y[i] = x[i]            for i constrained   (identity rows)
  residual  : r[i] = x[i] - g[i]     for i constrained
  assemble  : rows i replaced by identity, columns untouched
Applied once per sweep over the index vector; no branching in the inner
loop, no masking in the solver. zero_dirichlet_trace, is_dirichlet and
the #105 question "where do identity rows belong" disappear.

### 5.4 Kind and Args

The loop forwards both unchanged. Kind selects the local_operator
overload at compile time; Args is passed to apply, residual, update and
the reduce functor. The current loops map as:

    Elliptic::trace_to_flux        apply(Stiffness, {t})
    Elliptic::residual_flux        residual(Stiffness, {t})
    Elliptic::trace_to_flux_mat    assemble(Stiffness, {t})
    Parabolic::set_data            update({t})
    MassApproxEigenvalue           apply(Mass, {})
    ShiftedInverseEigenvalue       apply/assemble(Shifted, {sigma})
    Hyperbolic (ne9)               apply/residual(Stage, {t, k}), update
    Generic (ne9)                  the loop itself
    NonlinearEigenvalue, Newton    drivers on top of apply/residual

### 5.5 reduce

    loop.reduce(Kind, args, x, f, combine)

gathers the local trace, looks up the operator, calls
f(edge, op, lambda, args), folds with combine (default: sum, elementwise
for arrays). Under MPI sum/max are an allreduce; append maps to the
collective HDF5 write.

Examples:

    // Timoshenko dual error: reads State and analytic solution
    double err = loop.reduce(Stage{}, {t,0}, lam,
      [&](const Edge& e, const auto&, const auto& l, const auto& a)
      { return solver.dual_error_sq(e, l, a.time); });

    // Timoshenko energy: array of components, summed elementwise
    auto en = loop.reduce(Stage{}, {t,0}, lam,
      [&](const Edge& e, const auto&, const auto& l, const auto& a)
      { return solver.energy(e, l, a.time); });

    // Diffusion L2 error: needs the local solution from the operator
    double err = loop.reduce(Stiffness{}, {t}, lam,
      [&](const Edge& e, const auto& op, const auto& l, const auto& a)
      { return solver.l2_error_sq(e, op.local_solution(e, l, a), a.time); });

    // adaptive dt: largest local wave speed, combine = max
    // plotting / in-situ histograms: per-edge record, combine = append

Reductions differ in combine, in whether they need the operator, and in
result type. The loop knows none of this; the solver has no reduce member.

## 6. Parallelism

Shared memory
- apply/residual/reduce: parallel over edges; scatter by node coloring,
  per-thread buffers, or atomics (internal to the loop).
- update: parallel over edges, each edge owns its State.
- prepare: parallel over slots (§5.2).
- constraints: one pass over an index vector after the scatter.
- solver requirement: apply/residual are const, no lazy caches.

Distributed memory (distribute_domain exists for hyEdge_dim == 1)
- Vec is a distributed trace vector with ghost exchange on hypernodes.
- operators and store are rank-local; assemble yields rank-local triplets
  with global indices; ConstraintSet holds rank-local indices.
- reduce combines with an allreduce.
- LocalSolver and LocalOperator never see ranks.

Hybrid: the two layers compose because the local interface has no hidden
state.

## 7. What goes away

- five global loop classes, prototype_* macros, HAS_MEMBER_FUNCTION,
- HYPERHDG_LOCAL_TRACE_TO_FLUX_MAT and the hand-written
  Diffusion::trace_to_flux_mat (both become apply<k>),
- factorize_once and its sampling check (Stored + geometry key),
- lazy factorization in solve_local_problem (moves to local_operator),
- is_dirichlet, zero_dirichlet_trace, constrained() (ConstraintSet),
- matrix() (never introduced),
- reduce as a solver member (never introduced),
- the rule "the condensed matrix must not depend on hyper_edge.data"
  (true by construction: the operator holds no State).

## 8. Where encapsulation helps and where it hurts

Keep the wall:
- loop ↔ solver: one loop for every physics; parallel scheduling never
  learns the physics.
- solver ↔ cache lifetime: how/what-depends-on-what (solver) vs
  when/for-whom (loop); three code paths become one policy.
- State ↔ operator: a correctness tool; makes operators shareable and
  apply thread-safe.
- Kind/Args ↔ loop: new operators need no loop code.
- constraints ↔ solver: runtime policy, one index pass, no branching
  inside the operator.

Do not build the wall:
- inside the solver: no per-block apply objects; blocks are meaningless
  apart and have one client.
- geometry: pass the edge, read normals/integrals; copy nothing in.
- index mapping edge dofs ↔ global dofs: must stay visible to the loop.
- the operator: no abstract base, no virtual apply; `requires` + static
  dispatch, zero cost in a hot loop of tiny dense ops.
- the store: prepare/release and a map, not a policy class or a cache
  framework.

Rule: a wall where two sides evolve independently; none where only one
client exists and the pieces are meaningless apart.

## 9. Migration

1. Add LocalOperator + GlobalLoop beside Generic in ne9 with an adapter
   that wraps any existing solver's trace_to_flux/residual_flux as an
   Ephemeral operator and its is_dirichlet as a ConstraintSet predicate.
2. Port TimoshenkoWave: assemble_schur → local_operator, State keeps
   coefficients and trace only, remove zero_dirichlet_trace.
3. Port Diffusion with block apply (hoisted integrals, multi-rhs dgetrs)
   and a geometry key for UnitCube.
4. Python: one nanobind template over Kind/Args; sparse_stiff_mat etc.
   stay as thin aliases.
5. Retire the old loops.

## 10. Repository structure (planned, separate branch `structure` off main)

Independent of the loop redesign; not stacked on it. The only shared piece, NumPy across the
binding boundary (the ndarray overloads of Elliptic's methods in bind_python.hxx), is cut into
its own small branch `numpy-bindings` off main that both `operator-loop` and `structure` merge.

Layout:

    tests/
      unit/        qr_handedness.cxx, point_functions.cxx, dense_la.cxx   fast, no solver
      params/      manufactured solutions (parameter structs), shared by everything below
      solvers/
        <solver>.hxx   instantiation ladder: (dim, degree) types and names, one per solver
        <solver>.cxx   C++ smoke test and example: first ladder entry, coarse mesh, CG
                       (tests/solvers/cg.hxx, matrix-free, ~20 lines), error below a threshold
        <solver>.py    python ladder test: every entry, convergence rates
        bind.cxx       one generic nanobind module source, compiled once per solver with the
                       solver's ladder header injected; class names <Solver>_D<dim>_P<deg>
    examples/        python only: the just-in-time path (hyperhdg.load) on geo files, a few
                     scripts (diffusion, timoshenko); the C++ examples are replaced by
                     tests/solvers/*.cxx, which are compiled and run

CMake: foreach(solver): executable from solvers/<solver>.cxx, nanobind module from bind.cxx,
ctest running <solver>.py against that module (skipped without NumPy).

Rules:
- The ladder header is the single source of instantiations, read by the C++ test and the
  binding. Python knows only the class-name convention and writes its ladder as plain numbers;
  a mismatch fails with an attribute error naming the missing class. Nothing else crosses.
- Tests never compile: they import the prebuilt ladder modules. The JIT path lives in
  examples/ only (note at the top of tests/README.md).
- Ladders stay small (dims 1-2, degrees 1-3); 3D high degree is research, not test suite.
- Each <solver>.py is solve(opt) + a loop over the ladder, so the file reads as the example of
  the bindings (preallocate y, rows, cols, vals via n_triplets; no helpers hiding it) and runs
  as the test. The examples/ scripts share the same solve after the class is obtained.
- Convergence asserts: observed rate on the last step >= degree + 1 - 0.3, and the finest
  error below a bound. No equality on rates.
- Removed: tests_c++, tests_python (incl. the cython local_solvers/run.py), reproducibles_python
  (research scripts move to python/ next to the drivers), cython/, import/, examples/*.cxx.

Order: numpy-bindings (small, first) -> structure (on the old loops) -> after operator-loop
lands, a follow-up swaps the loop into tests/solvers and deletes python/operator_loop*.py.
