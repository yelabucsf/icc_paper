#!/bin/bash

out_dir="/data/analysis_outputs/QTL_analysis/chromHMM_enrichment"

# for cell_type in B_Memory B_Naive CD4T_Memory CD4T_Naive CD8T_Memory CD8T_Naive MAIT NK_dim Treg cDC2 cM gdT ncM;
# do
#     echo $cell_type
#     python chromHMM_caQTL_vs_eQTL_enrichment.py --celltype $cell_type -o $out_dir/${cell_type}_chromHMM_caQTL_vs_eQTL_enrichment
#     python chromHMM_activating_vs_repressive_enrichment.py --celltype $cell_type -o $out_dir/${cell_type}_chromHMM_activating_vs_repressive_enrichment
# done

for cell_type in B_Memory B_Naive CD4T_Memory CD4T_Naive CD8T_Memory CD8T_Naive MAIT NK_dim Treg cDC2 cM gdT ncM;
do
    echo $cell_type
    python chromHMM_caQTL_vs_eQTL_enrichment.py --celltype $cell_type -o $out_dir/${cell_type}_chromHMM_caQTL_vs_eQTL
done