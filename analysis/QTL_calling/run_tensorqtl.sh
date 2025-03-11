#!/bin/bash


# conda activate tensorQTL
mateqtl_dir="/data_vol/icc_objects/QTL_calling/rna"
cell_type="B_cell"
chr=1

for group in all_inds EUR EAS AFR AMR;
do
    mkdir ${mateqtl_dir}/outputs/${group}
    for num_pheno_pcs in 3 5 10;
    do
        for num_geno_pcs in 3 5 10;
        do
            python run_tensorqtl.py -g /data_vol/icc_objects/genotypes/vcfs/plink/chr${chr}.dose -p ${mateqtl_dir}/inputs/${group}/${group}_${cell_type}.phenos.tensorqtl.sorted.bed -c ${mateqtl_dir}/inputs/${group}/${group}_${cell_type}_${num_pheno_pcs}_phenoPCs_${num_geno_pcs}_genoPCs.cov -n ${mateqtl_dir}/outputs/${group}/${group}_${cell_type}_chr${chr}_${num_pheno_pcs}_phenoPCs_${num_geno_pcs}_genoPCs_tensorqtl --chrom ${chr}
        done
    done
done

