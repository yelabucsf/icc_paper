#!/bin/bash

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate aws

cells=( CD8T_Naive cDC1 cM MAIT ncM NK_dim Plasma Prolif B_Memory CD4T_Memory CD8T_Memory CD8T_platelet cDC2 gdT NK_bright pDC Progen Treg )


level='L2_manual'

mkdir -p /data2/processed/QTL_calling/atac/permissive/${level}/${cell}
echo $level

for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/  --exclude "*" --include "*cis_qtl_snp*.txt" ; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; echo 'annotate files' ; conda activate tensorQTL; python add_qtl_annots_atac.py -l ${level} -c ${cell} ; echo 'copy files to s3'; conda activate aws ; aws s3 sync /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/ s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ ; echo 'remove old files' ;  rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done


#cells2=['B_atypical','B_Naive','CD4T_Naive','CD8T_Naive','cDC1','cM','MAIT','ncM','NK_dim','Plasma','Prolif',
#        'B_Memory','CD4T_Memory','CD8T_Memory','CD8T_platelet','cDC2','gdT','MK_2','NK_bright','pDC','Progen','Treg']
#cells completed 'B_Naive','B_Memory''CD4T_Naive'

# copy files from s3


# run annotation script



#copy outputs to s3


