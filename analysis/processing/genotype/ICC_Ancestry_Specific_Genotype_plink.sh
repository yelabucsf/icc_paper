

plinkExec="/data/Genotype/plink2"
genotypedir="/data/Genotype/"

ancestries=("AFR" "AMR" "EAS" "EUR")

# function processGenotype {
#   for thisanc in ${ancestries[@]};
#     do
#         echo ${thisanc}
#         mkdir -p ${genotypedir}/Ancestry_Specific/${thisanc}
#         $plinkExec --pfile ${genotypedir}/chr${1}.dose \
#             --no-sex --allow-no-sex \
#             --keep ${genotypedir}/Ancestry_Specific/ICC_${thisanc}_individuals.txt --make-pgen --out ${genotypedir}/Ancestry_Specific/${thisanc}/chr${1}
#     done
# }
# export -f processGenotype
# parallel processGenotype ::: {1..22}


for thisChr in {1..22}
do
    for thisanc in ${ancestries[@]};
    do
        echo ${thisanc}
        mkdir -p ${genotypedir}/Ancestry_Specific/${thisanc}
        $plinkExec --vcf ${genotypedir}/vcf/chr${thisChr}.dose.vcf.gz \
            --no-sex --allow-no-sex \
            --keep ${genotypedir}/Ancestry_Specific/ICC_${thisanc}_individuals.txt --make-pgen --out ${genotypedir}/Ancestry_Specific/${thisanc}/chr${thisChr}
    done
done

### Parallel for each chromosome
function processVCF {
  for thisanc in ${ancestries[@]};
    do
        echo ${thisanc}
        mkdir -p ${genotypedir}/Ancestry_Specific/${thisanc}
        $plinkExec --vcf ${genotypedir}/vcf/chr${1}.dose.vcf.gz \
            --no-sex --allow-no-sex \
            --keep ${genotypedir}/Ancestry_Specific/ICC_${thisanc}_individuals.txt --make-pgen --out ${genotypedir}/Ancestry_Specific/${thisanc}/chr${1}
    done
}
export -f processVCF
parallel processVCF ::: {1..22}
