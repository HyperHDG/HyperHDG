#!/bin/bash
# Layouts of the local operator of the Timoshenko beam solver (BeamStiffness, the default, against
# Condensed and Factored) on a random mikado network: time of prepare / apply / residual /
# assemble, size of the store, agreement of the results. Run from the repository root:
#   experiments/oploop-01-timoshenko-layouts.sh [MASS] [R]
# The network comes from make_geo2.py (MAKE_GEO, by default that of the ne9 worktree); it needs
# HDF5 (.geo.h5). Output: output/oploop-01/.
set -euo pipefail
MASS=${1:-400}
R=${2:-0.05}
: ${MAKE_GEO:=$HOME/phd/HyperHDG/ne9/experiments/make_geo2.py}
: ${CXX:=g++}
OUT=output/oploop-01
mkdir -p $OUT
cp $0 experiments/oploop-01-timoshenko-layouts.cxx $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

DOM=$OUT/mikado.geo.h5
python3 $MAKE_GEO --mikado $MASS $R --seed 1 --dirichlet xmin=1 xmax=1 ymin=1 ymax=1 -o $DOM \
  | tee $OUT/make_geo.log

FLAGS="-std=c++20 -I. -Iinclude -Isubmodules/tensor_product_polynomials.git/include
  -Isubmodules/tensor_product_chain_complex.git/include -DHYPERHDG_HAVE_HDF5"
LIBS="-llapack -lhdf5"
: > $OUT/results.yaml
for DEG in 1 2 3; do
  $CXX $FLAGS -O3 -DNDEBUG -DDEGREE=$DEG experiments/oploop-01-timoshenko-layouts.cxx \
    -o $OUT/layouts_p$DEG $LIBS
  $OUT/layouts_p$DEG $DOM 20 | tee -a $OUT/results.yaml
done
# Once with assertions (rigid motions without flux, no coupling of shear and curvature, on the
# stiffnesses of the network), degree 1.
$CXX $FLAGS -O2 -DDEGREE=1 experiments/oploop-01-timoshenko-layouts.cxx -o $OUT/layouts_debug $LIBS
$OUT/layouts_debug $DOM 1 | tee $OUT/debug.yaml
