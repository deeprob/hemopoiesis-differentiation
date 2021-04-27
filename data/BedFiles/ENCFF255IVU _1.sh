#!/bin/bash

wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

bigBedToBed  ${1} ${1%%.bigBed}".bed"


cat ${1%%.bigBed}".bed" | awk 'BEGIN{OFS="\t"}{if (NR>1 && prev_chr==$1 && prev_chr_e<=$2) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > ${1%%.bigBed}"_n_o.bed"


cut -f1,2,3,4 ${1%%.bigBed}"_n_o.bed" | sort -k1,1 -k2,2n > ${1%%.bigBed}"_n_o_sort.bedGraph"


bedGraphToBigWig ${1%%.bigBed}"_n_o_sort.bedGraph" mm10.chrom.sizes ${1%%.bigBed}".bw"

cat ${1%%.bigBed}".bed" | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > ${1%%.bigBed}"_master_peak_list.withpkname.bed"


bigWigAverageOverBed ${1%%.bigBed}".bw" ${1%%.bigBed}"_master_peak_list.withpkname.bed" ${1%%.bigBed}"_master_peak_list.bigbedfile1.tab"


