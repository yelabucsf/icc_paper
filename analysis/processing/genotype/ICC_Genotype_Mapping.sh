
#data = fread("/data/person_demo_anc.2025-03-14_corrected.txt")
#write.table(data$V1, file = "/data/Genotype/ICC_individuals.txt", col.names = F, row.names = F, sep = "\n", quote = F)
# for(this.anc in c("AFR", "AMR", "EAS", "EUR")){
#     this.inds = data$V1[which(data$Assignment == this.anc)]
#     write.table(this.inds, file = paste0("/data/Genotype/ICC_individuals_", this.anc, ".txt"), col.names = F, row.names = F, sep = "\n", quote = F)
# }

plinkExec="/data/Genotype/plink2"
bcftoolsExc="/data/bcftools/bcftools"
genotypedir="/data/Genotype/"

ancestries=("AFR" "AMR" "EAS" "EUR")


for thisChr in {1..22}
do
    $bcftoolsExc reheader -s /data/vcf_mapping.2025-03-14.txt -o ${genotypedir}/vcf_reheader/chr${thisChr}_reheader.vcf.gz  ${genotypedir}/vcf/chr${thisChr}.dose.vcf.gz
    $plinkExec --vcf ${genotypedir}/vcf_reheader/chr${thisChr}_reheader.vcf.gz \
        --no-sex --allow-no-sex \
        --keep ${genotypedir}/ICC_individuals.txt --make-pgen --out ${genotypedir}/plink/chr${thisChr}
    for thisanc in ${ancestries[@]};
    do
        #mkdir -p ${genotypedir}/Ancestry_Specific/${thisanc}
        $plinkExec --vcf ${genotypedir}/vcf_reheader/chr${thisChr}_reheader.vcf.gz \
            --no-sex --allow-no-sex \
            --keep ${genotypedir}/ICC_individuals_${thisanc}.txt --make-pgen --out ${genotypedir}/Ancestry_Specific/${thisanc}/chr${thisChr}
    done
done

### Merged genotype file 
echo '/data/Genotype/vcf_reheader/chr1_reheader.vcf.gz' > ${genotypedir}/vcf_reheader_merge.txt
for thisChr in {2..22}
do
    echo '/data/Genotype/vcf_reheader/chr'${thisChr}'_reheader.vcf.gz' >> ${genotypedir}/vcf_reheader_merge.txt
done

$bcftoolsExc concat -f ${genotypedir}/vcf_reheader_merge.txt -Oz -o ${genotypedir}/vcf_reheader/combined_reheader.vcf.gz

$plinkExec --vcf ${genotypedir}/vcf_reheader/combined_reheader.vcf.gz \
    --no-sex --allow-no-sex \
    --keep ${genotypedir}/ICC_individuals.txt --make-pgen --out ${genotypedir}/plink/combined


###### Remove individuals with incorrect Redcap / Genotype mapping
sed -E -e "/CUPP_140|CUPP_185/d" ICC_individuals.txt > Cleaned/ICC_individuals_cleaned.txt

for thisanc in ${ancestries[@]};
    do
    sed -E -e "/CUPP_140|CUPP_185/d" ICC_individuals_${thisanc}.txt > Cleaned/ICC_individuals_${thisanc}_cleaned.txt
done

$plinkExec --pfile ${genotypedir}/plink/combined \
    --no-sex --allow-no-sex \
    --keep ${genotypedir}/Cleaned/ICC_individuals_cleaned.txt --make-pgen --out ${genotypedir}/plink/combined_cleaned

for thisChr in {1..22}
do
    $plinkExec --pfile ${genotypedir}/plink/chr${thisChr} \
        --no-sex --allow-no-sex \
        --keep ${genotypedir}/Cleaned/ICC_individuals_cleaned.txt --make-pgen --out ${genotypedir}/Cleaned/chr${thisChr}_cleaned
    for thisanc in ${ancestries[@]};
    do
        mkdir -p ${genotypedir}/Cleaned/Ancestry_Specific/${thisanc}
        $plinkExec --pfile ${genotypedir}/plink/chr${thisChr} \
            --no-sex --allow-no-sex \
            --keep ${genotypedir}/Cleaned/ICC_individuals_${thisanc}_cleaned.txt --make-pgen --out ${genotypedir}/Cleaned/Ancestry_Specific/${thisanc}/chr${thisChr}_cleaned
    done
done

