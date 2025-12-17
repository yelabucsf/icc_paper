### Script to LD clump eQTL hits for each eGene in each tissue

level=$1
tissue=$2 #loop over tissues

supp_data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/gwas_eqtl" # see Supplementary Data
lddir="/data_vol/FINAL/Genotypes/Cleaned"  # in-cohort LD
data_dir="/data_vol/FINAL/QTL/rna/${level}"
codedir="/data/icc_paper/analysis/QTL_analysis/qtl_properties_mostafavi/icc_eqtl_process"
outdirs="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/rna"
outdir=$outdirs/$level/$tissue
mkdir -p $outdir

#Get eGenes==========

outfile=$outdir/$tissue".eGenes.txt"
touch $outfile
for chr in {1..22}
do
    infile=$data_dir/$tissue/$tissue".cis_qtl_egene.${chr}.txt"    
    cat $infile | awk 'BEGIN{FS="\t"}($18<=0.05){print $1}' | sort | uniq >> $outfile
done

#format hg38 files to have correct ref and alt alleles=========
infile=$data_dir/$tissue/$tissue".cis_qtl_snp.all_significant.annot.txt"
outfile=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
cat $infile | awk 'BEGIN{FS="\t"}{print "chr"$1"_"$2"_"$4"_"$5"_b38",$7,$13,$3}' | tail -n +2 | awk '{print $1,$2,$3,$4}' > $outfile

infile=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
outfile=$outdir/$tissue".pos_hg38.bed"
cat $infile | awk 'BEGIN{FS="_"}{print $1,$2,$5}' | awk '{print $1,$2-1,$2,$6}' > $outfile

changefile=$supp_data_dir/"auxiliary_files/hg38ToHg19.over.chain.gz" 
maffile=$supp_data_dir/"snp_annotations/maf01.kgp.snps"
bimfile="/data_vol/FINAL/Genotypes/Cleaned/snps.txt"

infile=$outdir/$tissue".pos_hg38.bed"
outfile1=$outdir/$tissue".pos_hg19.bed"
outfile2=$outdir/$tissue".pos_hg38.failed"
liftover="/data/tools/liftOver"
$liftover $infile $changefile $outfile1 $outfile2 #using LiftOver software

infile1=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
infile2=$outdir/$tissue".pos_hg38.bed"
infile3=$outdir/$tissue".pos_hg19.bed"

script=$codedir/"keep_correct_alleles.R" #restrict to eQTLs that match KGP SNPs with MAF>0.01
Rscript $script $maffile $bimfile $infile1 $infile2 $infile3

#LD clumping==========

clumpdir=$outdir/"clump"
mkdir $clumpdir

cd $clumpdir
rm *

infile=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
outfile=$clumpdir/$tissue".signif_variant_gene_pairs.hg38.txt"
cp $infile $outfile
awk -F' ' '{print > $2".txt"}' $outfile

rm $outfile
cd ../

genelist=$outdir/$tissue".eGenes.txt"
script=$codedir/"get_clump_assoc.R"

concatfile=$outdir/$tissue".clumped.snps"
rm $concatfile

#loop over eGenes
while IFS= read -r line 
do

gene="`echo $line`"

infile=$clumpdir/$gene".txt"
chr="`head -n 1 $infile | awk '{print $1}'`"
tempfile=$clumpdir/$gene".temp.txt"
echo "CHR GENE SNP P" > $tempfile
cat $infile >> $tempfile
mv $tempfile $infile

ldfile=$lddir/"chr"$chr"_cleaned"
gwasfile=$clumpdir/$gene".txt"
outfile=$clumpdir/$gene
plink="/data/tools/plink_files/plink"
$plink --bfile $ldfile --clump $gwasfile --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out $outfile #using plink software

infile=$clumpdir/$gene".clumped"
temp=$clumpdir/$gene".clumped.temp"
Rscript $script $infile $temp
awk -v tissue=$tissue -v gene=$gene '{print tissue,gene,$0,NR}' $temp >> $concatfile
rm $temp

done < $genelist

maffile=$supp_data_dir/"snp_annotations/maf01.kgp.snps"
infile1=$outdir/$tissue".clumped.snps"
infile2=$outdir/$tissue".pos_hg38.bed"
infile3=$outdir/$tissue".pos_hg19.bed"
outfile=$outdir/$tissue".clumped.hg19.snps"

script=$codedir/"lift_snps_new.R" #lift snps to hg19
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile