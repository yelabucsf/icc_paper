#!/bin/bash
# Subset vcf genotypes to 333 samples in anndata

dir_dat="./dat/processed/cm_il1b/reqtl/vcfs"
mkdir -p $dir_dat

samples="./dat/processed/cm_il1b/reqtl/sample_ids.txt"

vcf_list=()

for chrm in $(seq 1 22) X; do
    echo "Chromosome $chrm"
    vcf_in="./dat/processed/imputed/oct_2024/chr${chrm}.dose.rehdr.uniq.r2flt0.3.bcf"
    vcf_out="${dir_dat}/chr${chrm}.dose.flt.bcf"

    bcftools view \
        -S "$samples" \
        -O u \
        "$vcf_in" | \
    bcftools +fill-tags \
        - \
        -- -t MAF,AF,AC,AN | \
    bcftools filter \
        -i 'INFO/MAF>=0.05' \
        -O b \
        -o "$vcf_out" \
        --threads 4
    bcftools index "$vcf_out"
    vcf_list+=("$vcf_out")
done

bcftools concat \
    -O b \
    -o "${dir_dat}/chrALL.dose.flt.bcf" \
    --threads 4 \
    "${vcf_list[@]}"
bcftools index "${dir_dat}/chrALL.dose.flt.bcf"
