#!/usr/bin/env bash
# 5-minute demo of the explicit nanobind binding path (see python/README.md).
# Run from the repo root: bash python/demo.sh
# Non-destructive: the live-edit step restores python/timowave_py.cxx on exit.

set -euo pipefail

step() { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
pause() { if [ -t 0 ]; then read -rp "-- enter to continue --"; fi; }

MODULE_CXX=python/timowave_py.cxx
SO=$(echo build/openblas/python/timowave_py.*.so)

restore() { git checkout -- "$MODULE_CXX"; }
trap restore EXIT

step "1. precompiled module in action: scipy driver, same YAML keys as the PETSc driver"
python python/timowave.py -deg 1 -nx 4 -nt 32 -T 1 -domain domains/single1.geo
pause

step "2. this is EVERYTHING that gets compiled (whole module source):"
cat "$MODULE_CXX"
pause

step "3. no petsc anywhere near it:"
ldd "$SO" | grep -i petsc || echo "  (no libpetsc in $SO)"
pause

step "4. marginal cost of a new instantiation: one line + incremental rebuild"
sed -i 's|.*TimoWave4_P2S2.*|&\n  HyperHDG::bind_python<HDGTimoWave<3, 2, Complex>>(m, "TimoWave4_P3S2");|' "$MODULE_CXX"
git diff --stat "$MODULE_CXX"
time cmake --build --preset openblas --target timowave_py >/dev/null
python - <<'EOF'
import sys; sys.path.insert(0, "build/openblas/python")
import timowave_py as hy
hdg = hy.TimoWave4_P3S2("domains/single1.geo", [1.0, 0.03125, 0.0])
print("TimoWave4_P3S2 constructed, size_of_system =", hdg.size_of_system())
EOF
pause

step "5. (live, in an editor) error experience"
cat <<'EOF'
Introduce a typo in a template argument in both systems and compare:
  - timowave_py.cxx: clangd flags it as you type; the build fails once, at build time.
  - legacy config string: accepted silently, fails at import time with raw compiler
    stderr in the middle of the script/notebook run.
EOF

step "done (timowave_py.cxx restored, module left with the extra instantiation until next build)"
