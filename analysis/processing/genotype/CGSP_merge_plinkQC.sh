# interact -n 4 -c 1 -p shared-name --t 6:00:00


cd /data/abattle4/taiboli/ICC/Genotype/CGSPLL240202

plinkExc="/home/tli64/data-abattle4/taiboli/plink"
kingExc="/home/tli64/data-abattle4/taiboli/king"


$plinkExc --merge-list mergeList.txt --make-bed --out Merged/CGSP_GDA01_MEGA2_5

###

$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5 --missing --out Merged/CGSP_GDA01_MEGA2_5
# 
# Delete SNPs with missingness >0.20 fist .
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5 --geno 0.2 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-removeSNPs_0.2
# Total genotyping rate is 0.901857.
# 326476 variants removed due to missing genotype data (--geno).
# 1661576 variants and 593 people pass filters and QC.


# Delete individuals with missingness >0.2.
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPs_0.2 --mind 0.2 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-removeSNPs_0.2-removeInds_0.2
# 5 people removed due to missing genotype data (--mind).
# Total genotyping rate in remaining samples is 0.997248.
# 1661576 variants and 588 people pass filters and QC.


###### Iteration 1: Use the same criteria as before

# Delete SNPs with missingness >0.02 fist .
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPs_0.2-removeInds_0.2 --geno 0.02 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-removeSNPs
# Total genotyping rate is 0.997248.
# 14153 variants removed due to missing genotype data (--geno).
# 1647423 variants and 588 people pass filters and QC.


# Delete individuals with missingness >0.05.
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPs --mind 0.05 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-removeSNPsInds
# 10 people removed due to missing genotype data (--mind).
# Total genotyping rate in remaining samples is 0.999177.
# 1647423 variants and 578 people pass filters and QC.



### Step 2: check sex
# plink --bfile plinkQC/MEGA2_5-removeSNPsInds --split-x hg38 no-fail --make-bed --out plinkQC/MEGA2_5-splitX
# ## data already contains an XY region, failed

$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPsInds --check-sex ycount --out Merged/CGSP_GDA01_MEGA2_5-checkSex
## analyze with R 

#### TODO 
## OK with sex check since all females have F < 0.31
# F statistic < 0, explanation https://www.animalgenome.org/community/angenmap/athread/13.php




### Step 3: MAF check
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPsInds --freq --out Merged/CGSP_GDA01_MEGA2_5-removeSNPsInds-MAFcheck

# Remove SNPs with a low MAF frequency.
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-removeSNPsInds --maf 0.01 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-postMAF
# Total genotyping rate is 0.999116.
# Total genotyping rate is 0.999177.
# 707184 variants removed due to minor allele threshold(s)
# 940239 variants and 578 people pass filters and QC.





### Step 4: hwe
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postMAF --hardy -out Merged/CGSP_GDA01_MEGA2_5-postMAF-HWE

$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postMAF --hwe 1e-4 --make-bed --out Merged/CGSP_GDA01_MEGA2_5-postHWE
# Total genotyping rate is 0.999106.
# --hwe: 76193 variants removed due to Hardy-Weinberg exact test.
# 864046 variants and 578 people pass filters and QC.





#
# ### Step 5: heterozygosity
# # heterozygosity of  admixed populations are higher https://www.biorxiv.org/content/10.1101/820241v1.full
#
# # prune for independent SNPs
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postHWE --indep-pairwise 50 10 0.1 --out Merged/CGSP_GDA01_MEGA2_5-postHWE-indep
# Pruning complete.  671192 of 862946 variants removed.




# plink --bfile plinkQC/MEGA2_5-postHWE --extract plinkQC/MEGA2_5-indep.prune.in --het --out plinkQC/MEGA2_5-indep-heterozygosity
#
# # try small sample approximation
# plink --bfile plinkQC/MEGA2_5-postHWE --extract plinkQC/MEGA2_5-indep.prune.in --het small-sample --out plinkQC/MEGA2_5-indep-heterozygosity-smallSample
# # 162287 variants and 175 people

#### Step 6: cryptic Relatedness
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postHWE --extract Merged/CGSP_GDA01_MEGA2_5-postHWE-indep.prune.in --genome --min 0.2 --out Merged/CGSP_GDA01_MEGA2_5-postHWE-pihat_min0.2
$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postHWE --extract Merged/CGSP_GDA01_MEGA2_5-postHWE-indep.prune.in --genome --out Merged/CGSP_GDA01_MEGA2_5-postHWE-genome-pihat



## http://people.virginia.edu/~wc9c/KING/manual.html
$kingExc -b Merged/CGSP_GDA01_MEGA2_5-postHWE.bed --related --prefix Merged/CGSP_GDA01_MEGA2_5-postHWE-king
$kingExc -b Merged/CGSP_GDA01_MEGA2_5-postHWE.bed --kinship --prefix Merged/CGSP_GDA01_MEGA2_5-postHWE-kinship ## this provided the kinship coefficient file
$kingExc -b Merged/CGSP_GDA01_MEGA2_5-postHWE.bed --ibdseg --prefix Merged/CGSP_GDA01_MEGA2_5-postHWE-ibdseg



####################

$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postHWE --missing --out Merged/CGSP_GDA01_MEGA2_5-postHWE


# cp ../../ICC/Genotype/CrossSectional/repeatFID-IIDs.txt ./
# echo -e "HCA-425\tHCA-425-2" >> repeatFID-IIDs.txt 
# echo -e "31_S3\tHCA_C09" 
# echo -e "65_S3\tHCA_C03"
# echo -e "66\tHCA_A02"

$plinkExc --bfile Merged/CGSP_GDA01_MEGA2_5-postHWE --remove CSGP_GDA_MEGA_removeRepeatRelated.txt --make-bed --out Merged/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd

# 864046 variants and 578 people pass filters and QC.

# removeRepeatRelated.txt


#####################################################
## Prepare input to imputation server
#####################################################

mkdir Merged/Impute/
cd Merged/Impute/

$plinkExc --bfile ../CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd --freq --out ../CGSP_GDA01_MEGA2_5-removeRepeatInd-MAFcheck

perl /home/tli64/data-abattle4/taiboli/ICC/ICC/Impute/HRC-1000G-check-bim.pl \
  -b ../CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd.bim \
  -f ../CGSP_GDA01_MEGA2_5-removeRepeatInd-MAFcheck.frq \
  -r /home/tli64/data-abattle4/taiboli/ICC/ICC/Impute/ALL.Variants.TOPMed.bravo.tab.gz -h \
  -v -n -o /home/tli64/data-abattle4/taiboli/ICC/Genotype/CGSPLL240202/Merged/Impute/Imputed

cd Imputed
sed 's/^/\/home\/tli64\/data-abattle4\/taiboli\//' Run-plink.sh > Run-plink-fixed.sh
bash Run-plink-fixed.sh
rm ../../TEMP*


ml parallel
function updateChr {
echo ${1}
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd-updated-chr"${1}".vcf > CGSP_GDA01_MEGA2_5-updated-chr"${1}"-withChr.vcf
}
export -f updateChr
parallel updateChr ::: {1..23}


ml htslib
function bgzipChr {
echo ${1}
bgzip -c CGSP_GDA01_MEGA2_5-updated-chr"${1}"-withChr.vcf > CGSP_GDA01_MEGA2_5-updated-chr"${1}"-withChr.vcf.gz
}
export -f bgzipChr
parallel bgzipChr ::: {1..23}






