#!/usr/bin/env sh
barcode="$1"
species="$2"
input=$(realpath $3)
log=$(realpath $4)
output=$(realpath $5)
container=$(realpath $6)
mkdir -p data/mykrobe/${barcode}
cd data/mykrobe/${barcode}
singularity exec ${container} mykrobe predict ${barcode} ${species} \
  --ont --seq ${input} --output ${output} 2> ${log}
