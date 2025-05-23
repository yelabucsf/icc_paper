#!/bin/bash

geno="/data_vol/icc_objects/genotypes/vcfs"
outdir="${geno}/plink"
mkdir -p $outdir

## create plink files
for chr in {1..22}; 
do 
    echo $i; 
    plink2 --vcf ${geno}/chr${chr}.dose.vcf.gz --out ${outdir}/chr${chr}.dose --threads 6 --set-all-var-ids @:# --make-pgen
done

chr="X"
plink2 --vcf ${geno}/chr${chr}.dose.vcf.gz --out ${outdir}/chr${chr}.dose --threads 6 --set-all-var-ids @:# --make-pgen
