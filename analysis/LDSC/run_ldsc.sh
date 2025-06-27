#!/bin/bash

## usage: ./run_ldsc.sh <ldsc_dir> <annot_list> <out_prefix> 

ldsc_dir=$1
annot_list=$2
out_prefix=$3
ldsc_ref_dir="/data_vol/References/ldsc_reference_files"

parallel -j 8 python2.7 /data/ldsc/ldsc.py --h2 ${ldsc_ref_dir}/sumstats/PASS_{}.sumstats --ref-ld-chr ${ldsc_ref_dir}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${annot_list} --w-ld-chr ${ldsc_ref_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --frqfile-chr ${ldsc_ref_dir}/1000G_Phase3_frq/1000G.EUR.QC. --out ${ldsc_dir}/out/{}_${out_prefix} ::: Alzheimer BMI1 Bipolar_Disorder Celiac Coronary_Artery_Disease Crohns_Disease Height1 IBD_deLange2017 Lupus Multiple_sclerosis Primary_biliary_cirrhosis Rheumatoid_Arthritis Schizophrenia Type_2_Diabetes Ulcerative_Colitis

