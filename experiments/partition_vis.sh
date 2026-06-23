#!/bin/bash
# partition_vis.sh — convert a morgan fiber network to .geo.h5, partition it
# with METIS, color the quotient graph with few colors, and render with net_vis.
#   make_geo.py (csv -> .geo.h5) -> partition_net.py -> net_vis.py
# All outputs land in output/partition_vis.<epoch>/ (symlinked as output/partition_vis).
set -xeo pipefail

: ${NET_DIR:=$HOME/phd/nextcloud/networks/morgan-2026-01-30/sca/net2}
: ${NPARTS:=256}
: ${COLOR_BY:=color}        # 'color' (quotient coloring) or 'partition' (raw index)
: ${PALETTE:=tab10}         # matplotlib colormap for the categorical colors
: ${FG:=black}              # color of category 0
: ${VIEW:=top}
: ${TUBE_RADIUS:=0}         # 0 = render as lines (fast, scale-independent)
: ${RES:=2000x2000}
: ${SHOW:=0}                # 1 = open interactive ParaView window
: ${OUTDIR:=output}
: ${PY:=venv/bin/python}
: ${PVPY:=pvpython}

NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
PNG=$OUT/$NAME.png

FN=submodules/fiber_network.git/fiber_network
MAKEGEO=$FN/make_geo.py
PARTITION=$FN/partition_net.py
NETVIS=$FN/net_vis.py
: ${METIS_LIB:=$(ls spack/openblas/.spack-env/._view/*/lib/libmetis.so 2>/dev/null | head -1)}

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cp $MAKEGEO $PARTITION $NETVIS $0 $OUT
git rev-parse HEAD > $OUT/ref
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/ref; fi

$PY $MAKEGEO -i $NET_DIR -o $DOMAIN --quirk morgan-2026-01-30
$PY $PARTITION $DOMAIN -n $NPARTS --metis-lib $METIS_LIB

# categories 1..N-1 (N = distinct values colored by); value 0 takes --fg / palette[0]
if [ "$COLOR_BY" = partition ]; then
  NCAT=$NPARTS
else
  NCAT=$($PY -c "import h5py; print(int(h5py.File('$DOMAIN','r')['domain'].attrs['n_colors']))")
fi
if [ "$NCAT" -gt 1 ]; then CATS=$(seq -s, 1 $((NCAT - 1))); else CATS=; fi

$PVPY $NETVIS $DOMAIN \
  --color-by $COLOR_BY --color-categories "$CATS" --color-palette $PALETTE --fg $FG \
  --warp-by none --arrows 0 \
  --view $VIEW -r $TUBE_RADIUS --resolution $RES -o $PNG --show $SHOW

set +x
echo "wrote $PNG"
if [ -n "${DISPLAY:-}" ] && command -v sxiv >/dev/null 2>&1; then sxiv "$PNG" & fi
