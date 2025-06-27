#!/bin/bash

plink2="/data/tools/plink2_files/plink2"
data_dir="/data_vol/FINAL/Genotypes/Cleaned"
for chr in {1..22}
do
    $plink2 --pgen $data_dir/chr${chr}_cleaned.pgen --pvar $data_dir/chr${chr}_cleaned.pvar --psam $data_dir/chr${chr}_cleaned.psam --make-bed --out $data_dir/chr${chr}_cleaned
    awk '{ $2 = $1":"$4":"$5":"$6; print }' OFS='\t' $data_dir/chr${chr}_cleaned.bim > $data_dir/chr${chr}_cleaned.fix_ids.bim
    mv $data_dir/chr${chr}_cleaned.fix_ids.bim $data_dir/chr${chr}_cleaned.bim
done