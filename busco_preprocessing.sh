#!/bin/bash
seqfile=$1
annfile=$2

gt gff3 -sort -tidy -retainids $annfile > tmp.gff3
gt extractfeat -type CDS -translate -join -retainids -seqfile $seqfile -matchdescstart < tmp.gff3
rm tmp.gff3 &> /dev/null
