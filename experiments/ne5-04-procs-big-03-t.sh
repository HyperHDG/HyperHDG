#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-04-procs-big}" # set default if unset
jq -sc 'map(.Stdout | {p: .mpi.sz, t_t2f}) | sort_by(.p).[]' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_ksp}) | sort_by(.p).[]' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_ksp}) | sort_by(.p).[]' $OUT.json
jq -sc 'map(.Stdout | {p: .mpi.sz, t_it: .t_iteration}) | sort_by(.p).[]' $OUT.json
