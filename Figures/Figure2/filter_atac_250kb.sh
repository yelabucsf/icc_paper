#!/bin/bash

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate aws


cells=( B_Naive CD4T_Naive CD8T_Naive cM MAIT ncM NK_dim B_Memory CD4T_Memory cDC2 gdT Treg )

cells=( CD4T_Naive CD8T_Naive cM MAIT ncM NK_dim B_Memory CD4T_Memory cDC2 gdT Treg )

level='L2_manual'
mkdir -p /data/icc_paper/Figures/Figure2/manhattan/${level}/${cell}
echo $level
for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; conda activate aws ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/  --exclude "*" --include "*.annot.txt" ; echo 'filter 250kb' ; for i in {1..22} ; do echo 'chr:'${i} ;  mkdir -p /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/ ; conda activate tensorQTL ; python /data/icc_paper/Figures/Figure2/filter_atac_outputs.py  -f /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/${cell}.cis_qtl_snp.${i}.annot.txt  -c 250000 -o /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/${cell}.cis_qtl_snp.${i}.annot.txt; cp /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/*egene* /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/; done ; conda activate aws ; aws s3 sync /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/ s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/cis_250kb/ ; echo 'remove old files' ; rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done




