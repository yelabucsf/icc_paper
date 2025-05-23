import genetics_plotting_utils as gu
import gc

#/data/processed/QTL_calling/atac/permissive/L1_manual/B/

cells2=['B_atypical','B_Naive','CD4T_Naive','CD8T_Naive','cDC1','cM','MAIT','ncM','NK_dim','Plasma','Prolif',
        'B_Memory','CD4T_Memory','CD8T_Memory','CD8T_platelet','cDC2','gdT','MK_2','NK_bright','pDC','Progen','Treg']


print('start analysis L2')
for cell in cells2:
    print(cell)
    for i in range(1, 23):
        print(i)
        qtl_snp_output_annot=gu.process_raw_tensorqtl_new('/data/processed/march2025/geneid_mapping.txt', str('/data/processed/genotypes/plink/chr'+str(i)+'_cleaned.pvar'), str('/data/processed/QTL_calling/atac/permissive/L2_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.txt'), str('/data2/processed/QTL_calling/atac/permissive/L2_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.annot.txt'), 1000000, str('/data/processed/genotypes/plink/chr'+str(i)+'_merged_freqs.afreq') )
        qtl_snp_output_annot=0
        gc.collect()
print('done')

##for cell in cells:
##    print(cell)
##    for i in range(22, 23):
##        print(i)
##        qtl_snp_output_annot=gu.process_raw_tensorqtl_dask('/data/processed/march2025/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots_filtered.h5ad', str('/data/processed/genotypes/plink/chr'+str(i)+'_cleaned.pvar'), str('/data/processed/QTL_calling/atac/permissive/L1_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.txt'), str('/data2/processed/QTL_calling/atac/permissive/L1_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.annot.txt'), [str('/data/processed/genotypes/plink/chr'+str(i)+'_ICC_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_AFR_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_AMR_freq.afreq'), str('/data/processed/genotypes/plink/chr'+str(i)+'_EAS_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_EUR_freq.afreq')], ['ICC_AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF'])
##
##cells2=['B_atypical','B_Naive','CD4T_Naive','CD8T_Naive','cDC1','cM','MAIT','ncM','NK_dim','Plasma','Prolif',
##        'B_Memory','CD4T_Memory','CD8T_Memory','CD8T_platelet','cDC2','gdT','MK_2','NK_bright','pDC','Progen','Treg']
##
##
##print('start analysis L2')
##for cell in cells2:
##    print(cell)
##    for i in range(1, 23):
##        print(i)
##        qtl_snp_output_annot=gu.process_raw_tensorqtl_dask('/data/processed/march2025/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots_filtered.h5ad', str('/data/processed/genotypes/plink/chr'+str(i)+'_cleaned.pvar'), str('/data/processed/QTL_calling/atac/permissive/L2_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.txt'), str('/data2/processed/QTL_calling/atac/permissive/L2_manual/'+str(cell)+'/'+str(cell)+'.cis_qtl_snp.'+str(i)+'.annot.txt'), [str('data/processed/genotypes/plink/chr'+str(i)+'_ICC_freq.afreq'), str('data/processed/genotypes/plink/chr'+str(i)+'_AFR_freq.afreq'), str('data/processed/genotypes/plink/chr'+str(i)+'_AMR_freq.afreq'), str('data/processed/genotypes/plink/chr'+str(i)+'_EAS_freq.afreq'), str('data/processed/genotypes/plink/chr'+str(i)+'_EUR_freq.afreq')], ['ICC_AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF'])
##
##        
##
