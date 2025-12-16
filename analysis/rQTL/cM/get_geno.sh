#!/bin/bash

dir_out="./dat/processed/imputed/final/"
mkdir -p "$dir_out"

aws s3 cp \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/Genotypes/combined.pgen \
    $dir_out

aws s3 cp \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/Genotypes/combined.pvar \
    $dir_out

aws s3 cp \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/Genotypes/combined.psam \
    $dir_out

dir_vcf="./dat/processed/imputed/oct_2024"
for chr in $(seq 1 22) X; do
    aws s3 cp \
        s3://gracie-icc/production_data/genotypes/processed/imputed/oct_2024/vcfs_uniq_r2flt0.3/chr${chr}.dose.rehdr.uniq.r2flt0.3.bcf \
        $dir_vcf
    aws s3 cp \
        s3://gracie-icc/production_data/genotypes/processed/imputed/oct_2024/vcfs_uniq_r2flt0.3/chr${chr}.dose.rehdr.uniq.r2flt0.3.bcf.csi \
        $dir_vcf
done
