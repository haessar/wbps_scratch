#!/bin/bash
input=$1
db=$2

input_fn=$(basename -- "$input")
db_fn=$(basename -- "$db")
omamer_fn="${input_fn%.*}.omamer"
output_dir="omark_${db_fn%.*}"

mkdir -p $output_dir
omamer search --db $db --query $input --out $output_dir/$omamer_fn
omark -f $output_dir/$omamer_fn -d $db -o $output_dir
