!#/bin/bash

conda activate tensorQTL

geno_dir='/data/processed/genotypes/plink/'
ancestry_files=( 'ICC_individuals_AFR_cleaned.txt',
	'ICC_individuals_AMR_cleaned.txt',
	'ICC_individuals_EAS_cleaned.txt',
	'ICC_individuals_EUR_cleaned.txt')
ancestry=( 'AFR' 'AMR' 'EAS' 'EUR' )
#chr14_cleaned


for chr in {1..22}; do echo $chr; plink2 --pfile ${geno_dir}chr${chr}_cleaned --freq --out ${geno_dir}chr${chr}_ICC_freq ; done


for chr in {1..22}; do echo $chr; for anc in ${ancestry[@]}; do echo $anc ; plink2 --pfile ${geno_dir}chr${chr}_cleaned --keep ${geno_dir}ICC_individuals_${anc}_cleaned.txt --freq --out ${geno_dir}chr${chr}_${anc}_freq ; done

