#!/bin/bash

## usage: ./make_annots <ldsc_dir> <bed_file> <out_file_prefix>

ldsc_dir=$1
bed_file=$2
out_file_prefix=$3

echo "LDSC dir: ${ldsc_dir}"
echo "Bed file: ${bed_file}"
echo "Out file prefix: ${out_file_prefix}"

echo "using hg19"
bimfile="/data_vol/References/ldsc_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim"
bfile="/data_vol/References/ldsc_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}"
snpfile="/data_vol/References/ldsc_reference_files/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.{}.snp"

# conda activate py27
liftOver ${bed_file}.bed /data/tools/hg38ToHg19.over.chain.gz ${bed_file}.hg19.bed ${bed_file}.unmap
parallel -j 8 python2.7 /data/ldsc/make_annot.py --bed-file ${bed_file}.hg19.bed --bimfile ${bimfile} --annot-file ${ldsc_dir}/annots/${out_file_prefix}.{}.annot.gz ::: {1..22}
parallel -j 8 python2.7 /data/ldsc/ldsc.py --l2 --bfile ${bfile} --ld-wind-cm 1 --annot ${ldsc_dir}/annots/${out_file_prefix}.{}.annot.gz --thin-annot --out ${ldsc_dir}/annots/${out_file_prefix}.{} --print-snps ${snpfile} ::: {1..22}
