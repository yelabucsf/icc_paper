#!/bin/bash

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate aws


#run pbmc
cells=( ICC )
level='Cohort'
#/data2/processed/QTL_calling/atac/permissive/Cohort/ICC/cis_250kb/
# test with 250kb windows

for cell in ${cells[@]} ; do echo $cell;  echo 'make plot' ; conda activate tensorQTL; echo 'plot rna';  python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot_rna.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; done

#for cell in ${cells[@]} ; do echo $cell;  echo 'make plot' ; conda activate tensorQTL; python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'make basic' ; python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'plot rna';  python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot_rna.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'make basic manhattan';  python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; done


##test rna only
#for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; conda activate aws ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/  --exclude "*" --include "*annot.txt" ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/ /data2/processed/QTL_calling/rna/${level}/${cell}/  --exclude "*" --include "*annot.txt"; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; for i in {1..22} ; do aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/rna/${level}/${cell}/ ; aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/ ; done; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; echo 'make plot' ; conda activate tensorQTL;  echo 'plot rna';  python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot_rna.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'make basic manhattan';  python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; done # echo 'remove old files' ; rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done


##real code
#for cell in ${cells[@]} ; do echo $cell; echo 'copy files to volume'; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/ /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/  --exclude "*" --include "*annot.txt" ; aws s3 sync s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/ /data2/processed/QTL_calling/rna/${level}/${cell}/  --exclude "*" --include "*annot.txt"; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; for i in {1..22} ; do aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/rna/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/rna/${level}/${cell}/ ; aws s3 cp s3://gracie-icc/production_data/analysis_outputs/FINAL/QTL/atac/permissive/${level}/${cell}/${cell}.cis_qtl_egene.${i}.txt /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/ ; done; ls /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/; echo 'make plot' ; conda activate tensorQTL; python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'plot rna';  python /data/icc_paper/Figures/Figure2/annotated_manhattan_plot_rna.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; echo 'make basic manhattan';  python /data/icc_paper/Figures/Figure2/basic_manhattan_plot.py -l ${level} -c ${cell} -d /data2/processed/QTL_calling/ ; done # echo 'remove old files' ; rm /data2/processed/QTL_calling/atac/permissive/${level}/${cell}/* ; done 
#



