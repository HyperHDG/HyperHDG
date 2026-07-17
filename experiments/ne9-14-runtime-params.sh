#!/bin/bash
# ne9-14: what do RUNTIME parameter functions cost?  wave4 (compile-time statics, inlined into
# the quadrature loops) vs wave4fn (identical formulas behind static inline std::function
# members -- one indirect call per quadrature point; see TestTimoWave4Fn in parameters.hxx).
# On the grid domain wave4 is not a manufactured solution (arm_axis assumes on-axis points), so
# e_* are meaningless here -- but they must AGREE between the two arms (identical fp ops), which
# is the correctness cross-check.  The numbers that matter:
#   t_t2f          matrix assembly (once; no parameter functions -> expect ~1)
#   t_rf / nt      residual assembly per step: RHS collocation integrals evaluate the parameter
#                  functions at every quadrature point -- the worst-case relative overhead,
#                  since wave4's functions are trivial
#   t_Timestepping / nt   full step incl. KSP solve -- the number the decision hangs on
# DECISION RULE (recorded up front): if the t_Timestepping overhead of wave4fn stays within a
# few percent, the future GlobalLoop::Generic configure() may take runtime parameter objects;
# if only t_rf moves while t_Timestepping is solve-dominated, likewise acceptable; otherwise
# parameters stay compile-time.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
LOGG=$OUT/make_geo.log
: ${NT:=20}
: ${REPS:=3}
NET="-pc_type net2as -net2as_p 7 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
KSP="-ksp_rtol 1e-9"
# serial on purpose: single-rank stage timings, no MPI noise
if [ -z "${PYTHON:-}" ]; then
  PYTHON=spack/$PRESET/.spack-env/view/bin/python
  $PYTHON -c 'import numpy, h5py' >/dev/null 2>&1 || PYTHON=python
fi
export PYTHON

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-19 grid convention: unit square, full clamp (type 63) on the boundary ring
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 1e-3"
for n in 65 129; do
  $PYTHON experiments/make_geo2.py --grid $n $DIR -o $DOMAIN-grid$n.geo.h5 | tee -a $LOGG
  $PYTHON experiments/gortz_constants.py $DOMAIN-grid$n.geo.h5 --mu --csv $OUT/gortz-grid$n.csv \
    | tee $OUT/gortz-grid$n.txt
done

# (deg, grid): deg 1 on the large grid, deg 5 (stages 1, real build) on the smaller one
run() { # $1 test  $2 deg  $3 grid  $4 rep
  local log=$OUT/$1-deg$2-grid$3-r$4.yaml
  $BUILD/timowave -test $1 -domain $DOMAIN-grid$3.geo.h5 -deg $2 -stages 1 \
    -nt $NT -T 1 $KSP $NET | tee $log
}
for rep in $(seq 1 $REPS); do
  for cfg in "1 129" "5 65"; do
    set -- $cfg
    run wave4   $1 $2 $rep
    run wave4fn $1 $2 $rep
  done
done

# summary table: per (deg, grid, rep) the three timings of both arms and their ratio
python - "$OUT" <<'EOF' | tee $OUT/summary.txt
import glob, re, sys, os
out = sys.argv[1]
def load(path):
    keys = {}
    for line in open(path):
        m = re.match(r'^(t_t2f|t_rf|t_Timestepping|e_rel|e_trace):\s*([0-9.eE+-]+)', line)
        if m: keys[m.group(1)] = float(m.group(2))
    return keys
runs = {}
for path in sorted(glob.glob(os.path.join(out, 'wave4*-deg*-grid*-r*.yaml'))):
    m = re.match(r'(wave4fn|wave4)-deg(\d+)-grid(\d+)-r(\d+)\.yaml', os.path.basename(path))
    runs[(m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)))] = load(path)
print(f"{'deg':>4} {'grid':>5} {'rep':>4} {'metric':>16} {'static':>12} {'fn':>12} {'fn/static':>10}")
for (test, deg, grid, rep), vals in sorted(runs.items()):
    if test != 'wave4': continue
    fn = runs.get(('wave4fn', deg, grid, rep), {})
    assert abs(vals.get('e_rel', 0) - fn.get('e_rel', 1)) < 1e-12, \
        f"e_rel mismatch deg{deg} grid{grid} r{rep}: identical fp ops expected"
    for key in ('t_t2f', 't_rf', 't_Timestepping'):
        a, b = vals.get(key), fn.get(key)
        if a and b:
            print(f"{deg:>4} {grid:>5} {rep:>4} {key:>16} {a:12.4e} {b:12.4e} {b/a:10.3f}")
EOF
