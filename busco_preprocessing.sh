#!/bin/bash
seqfile=$1
annfile=$2

# Clean seqids corrupted by BRAKER
sed 's/_length=[0-9]*//g' $annfile | gt gff3 -sort -tidy -retainids > tmp.gff3
gt extractfeat -type CDS -translate -join -retainids -seqfile $seqfile -matchdescstart < tmp.gff3
rm tmp.gff3 &> /dev/null
