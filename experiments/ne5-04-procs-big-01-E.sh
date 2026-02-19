#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-04-procs-big}" # set default if unset
jq -sc 'map(.Stdout | {p: .mpi.sz, t_t2f}) | sort_by(.p) | .[0].t_t2f as $t2f | .[] | {p, stage: "t2f", E: ($t2f / (.t_t2f * .p))}' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_ksp}) | sort_by(.p) | .[]' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_ksp}) | sort_by(.p) | .[0].t_ksp as $ksp | .[] | {p, stage: "ksp", E: ($ksp / (.t_ksp * .p))}' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_it: .t_iteration}) | sort_by(.p) | .[0].t_it as $it | .[] | {p, stage: "it", E: ($it / (.t_it * .p))}' $OUT.json

