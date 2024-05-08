#!/bin/bash
ctg_num=$1

samtools faidx ancylostoma_ceylanicum.PRJNA72583.WBPS19.genomic.fa ANCCEYDFT_Contig${ctg_num} > ANCCEYDFT_Contig${ctg_num}.fa
awk "/ANCCEYDFT_Contig$ctg_num\_/" braker.gff3 > ANCCEYDFT_Contig${ctg_num}_braker.gff
awk -v num="$ctg_num" '$1=="ANCCEYDFT_Contig"num' Ancylostoma_ceylanicum_helixer_hq.gff3 > ANCCEYDFT_Contig${ctg_num}_helixer.gff
awk -v num="$ctg_num" '$1=="ANCCEYDFT_Contig"num' ancylostoma_ceylanicum.PRJNA72583.WBPS19.annotations.gff3 > ANCCEYDFT_Contig${ctg_num}_wbps.gff

~/Artemis/art -Dblack_belt_mode=false ANCCEYDFT_Contig${ctg_num}.fa + ANCCEYDFT_Contig${ctg_num}_wbps.gff + ANCCEYDFT_Contig${ctg_num}_helixer.gff + ANCCEYDFT_Contig${ctg_num}_braker.gff