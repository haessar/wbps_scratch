#!/bin/bash
ctg_num=$1

samtools faidx schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa SM_V10_${ctg_num} > SM_V10_${ctg_num}.fa
awk "/SM_V10_$ctg_num\_/" Schistosoma_mansoni_braker3_full.gff3 > SM_V10_${ctg_num}_braker.gff
awk -v num="$ctg_num" '$1=="SM_V10_"num' Schistosoma_mansoni_helixer_full.gff3 > SM_V10_${ctg_num}_helixer.gff
awk -v num="$ctg_num" '$1=="SM_V10_"num' schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3 > SM_V10_${ctg_num}_wbps.gff

~/Artemis/art -Dblack_belt_mode=false SM_V10_${ctg_num}.fa + SM_V10_${ctg_num}_wbps.gff + SM_V10_${ctg_num}_braker.gff + SM_V10_${ctg_num}_helixer.gff
