#!/bin/bash

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate aws


cells=( B_Naive CD4T_Naive CD8T_Naive cM MAIT ncM NK_dim B_Memory CD4T_Memory cDC2 gdT Treg )
level='L2_manual'
cells=( CD8T_Naive cM Treg )

mkdir -p /data/icc_paper/Figures/Figure2/manhattan/${level}/${cell}
echo $level
for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; conda activate aws ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/cis_250kb/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/  --exclude "*" --include "*annot.txt" ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/ /data2/processed/QTL_calling/rna/${level}/${cell}/  --exclude "*" --include "*annot.txt"; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; for i in {1..22} ; do aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/rna/${level}/${cell}/ ; aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/ ;  cp /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/*egene* /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/; done ;echo 'make plot' ; conda activate tensorQTL; python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/  ; echo 'make basic manhattan';  python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ;  echo 'remove old files' ; rm /data2/processed/QTL_calling/rna/${level}/${cell}/* ; rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done


#for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; conda activate aws ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/  --exclude "*" --include "*annot.txt" ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/ /data2/processed/QTL_calling/rna/${level}/${cell}/  --exclude "*" --include "*annot.txt"; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; for i in {1..22} ; do aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/rna/${level}/${cell}/ ; aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/ ;  echo 'filter 250kb';  mkdir -p /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/ ; /data/icc_paper/Figures/Figure2/filter_atac_outputs.py  -f /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/${cell}.cis_qtl_snp.${i}.annot.txt  -c 250000 -o /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/${cell}.cis_qtl_snp.${i}.annot.txt; cp /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/*egene* /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/cis_250kb/; done ;echo 'make plot' ; conda activate tensorQTL; python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'plot rna';  python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot_rna.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'make basic manhattan';  python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ;  echo 'remove old files' ; rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done 

#cells2=['B_atypical','B_Naive','CD4T_Naive','CD8T_Naive','cDC1','cM','MAIT','ncM','NK_dim','Plasma','Prolif',
#        'B_Memory','CD4T_Memory','CD8T_Memory','CD8T_platelet','cDC2','gdT','MK_2','NK_bright','pDC','Progen','Treg']
#cells completed 'B_Naive','B_Memory''CD4T_Naive'

# copy files from s3


# run annotation script



#copy outputs to s3


