#!/usr/bin/env bash
# Same as run_mageck.sh, but on guides with median raw input count >= 100 (input_rep1-3),
# matching the paper's exclusion of sgRNAs with <100x coverage at the initial time point.
# Filtering happens BEFORE MAGeCK so dropped guides do not enter normalization or BH.
# Run from the repo root with MAGeCK 0.5.9.5 on PATH.
set -euo pipefail
O=mageck/og/ABE; K=$O/OG_ABE_COUNTS_input100.txt; F=$O/input100; C=input_rep1,input_rep2,input_rep3
declare -a TIS=(d15 spleen bonemarrow meninges)
declare -a TRT=(d15_rep1,d15_rep2,d15_rep3 spleen1,spleen2,spleen3,spleen4
                bonemarrow1,bonemarrow2,bonemarrow3,bonemarrow4,bonemarrow5
                meninges1,meninges2,meninges3,meninges4,meninges5)
mkdir -p $F/default $F/nc $F/logs
for i in "${!TIS[@]}"; do c=${TIS[$i]}; t=${TRT[$i]}
  mageck test -k $K -t $t -c $C -n $F/default/ABE_${c}_input > $F/logs/default_$c.log 2>&1 &
  mageck test -k $K -t $t -c $C --control-sgrna $O/control_sgrnas.txt --norm-method control \
    -n $F/nc/ABE_${c}_input > $F/logs/nc_$c.log 2>&1 &
done; wait
