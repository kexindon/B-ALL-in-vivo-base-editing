#!/usr/bin/env bash
# MAGeCK sgRNA tests for ABE OG, each tissue vs input_rep1-3 (same baseline as the empirical LFC-FDR table).
#   default  : median normalization
#   nc       : --control-sgrna + --norm-method control (negative-control normalization, as in s41587-026-03077-5)
#   ctrlnull : --control-sgrna only; sgRNA table is identical to default (controls only enter the gene-level RRA null)
# Run from the repo root with MAGeCK 0.5.9.5 on PATH.
set -euo pipefail
K=counts/og/OG_ABE_COUNTS.txt; O=mageck/og/ABE; C=input_rep1,input_rep2,input_rep3
declare -a TIS=(d15 spleen bonemarrow meninges)
declare -a TRT=(d15_rep1,d15_rep2,d15_rep3 spleen1,spleen2,spleen3,spleen4
                bonemarrow1,bonemarrow2,bonemarrow3,bonemarrow4,bonemarrow5
                meninges1,meninges2,meninges3,meninges4,meninges5)
mkdir -p $O/default $O/nc $O/ctrlnull $O/logs
for i in "${!TIS[@]}"; do c=${TIS[$i]}; t=${TRT[$i]}
  mageck test -k $K -t $t -c $C -n $O/default/ABE_${c}_input > $O/logs/default_$c.log 2>&1
  mageck test -k $K -t $t -c $C --control-sgrna $O/control_sgrnas.txt --norm-method control \
    -n $O/nc/ABE_${c}_input > $O/logs/nc_$c.log 2>&1
done
mageck test -k $K -t ${TRT[2]} -c $C --control-sgrna $O/control_sgrnas.txt \
  -n $O/ctrlnull/ABE_bonemarrow_input > $O/logs/ctrlnull.log 2>&1
