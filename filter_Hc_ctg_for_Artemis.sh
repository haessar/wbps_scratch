#!/bin/bash
ctg_num=$1

samtools faidx haemonchus_contortus.PRJEB506.WBPS19.genomic.fa hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon > hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon.fa
awk "/hcontortus\_chr$ctg_num\_/" Haemonchus_contortus_braker3_full.gff3 > hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_braker.gff
awk -v num="$ctg_num" '$1=="hcontortus_chr"num"_Celeg_TT_arrow_pilon"' Haemonchus_contortus_helixer_full.gff3 > hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_helixer.gff
awk -v num="$ctg_num" '$1=="hcontortus_chr"num"_Celeg_TT_arrow_pilon"' haemonchus_contortus_gca000469685v2.LT.fixed.gff3 > hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_anno.gff
awk -v num="$ctg_num" '$1=="hcontortus_chr"num"_Celeg_TT_arrow_pilon"' haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3 > hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_wbps.gff

~/Artemis/art -Dblack_belt_mode=false hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon.fa \
        + hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_wbps.gff \
        + hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_braker.gff \
        + hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_helixer.gff \
        + hcontortus_chr${ctg_num}_Celeg_TT_arrow_pilon_anno.gff
