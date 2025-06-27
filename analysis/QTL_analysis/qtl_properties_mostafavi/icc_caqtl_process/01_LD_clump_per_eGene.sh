### Script to LD clump eQTL hits for each eGene in each tissue

tissue=$1 #loop over tissues

supp_data_dir="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/gwas_eqtl" # see Supplementary Data
lddir="/data_vol/References/ldsc_reference_files/1000G_EUR_Phase3_plink" #1000G data used as LD ref panel because I don't have access to UKB
data_dir="/data_vol/FINAL/QTL/atac/permissive/L2_manual"
codedir="/data/icc_paper/analysis/QTL_analysis/qtl_properties_mostafavi/icc_caqtl_process"
outdirs="/data_vol/analysis_outputs/QTL_analysis/qtl_properties_mostafavi/gwas_eqtl/icc_caqtl_process"
outdir=$outdirs/$tissue
mkdir -p $outdir

#Get eGenes==========

outfile=$outdir/$tissue".eGenes.txt"
touch $outfile
for chr in {1..22}
do
    infile=$data_dir/$tissue/$tissue".cis_qtl_egene.${chr}.txt"    
    cat $infile | awk 'BEGIN{FS="\t"}($18<=0.05){print $1}' | sort | uniq >> $outfile
done

#lift hg38 stats file to hg19==========

infile=$data_dir/$tissue/$tissue".cis_qtl_pairs.significant_annotated.tsv"
outfile=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
cat $infile | awk 'BEGIN{FS="\t"}{print $12"_"$13"_"$14"_"$15"_b38",$1,$8,$2}' | tail -n +2 | awk '{print $1,$2,$3,$4}' > $outfile

infile=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
outfile=$outdir/$tissue".pos_hg38.bed"
cat $infile | awk 'BEGIN{FS="_"}{print $1,$2,$5}' | awk '{print $1,$2-1,$2,$6}' > $outfile

changefile=$supp_data_dir/"auxiliary_files/hg38ToHg19.over.chain.gz" 
maffile=$supp_data_dir/"snp_annotations/maf01.kgp.snps"

infile=$outdir/$tissue".pos_hg38.bed"
outfile1=$outdir/$tissue".pos_hg19.bed"
outfile2=$outdir/$tissue".pos_hg38.failed"
liftover="/data/tools/liftOver"
$liftover $infile $changefile $outfile1 $outfile2 #using LiftOver software

infile1=$outdir/$tissue".signif_variant_gene_pairs.hg38.txt"
infile2=$outdir/$tissue".pos_hg38.bed"
infile3=$outdir/$tissue".pos_hg19.bed"
outfile=$outdir/$tissue".signif_variant_gene_pairs.hg19.txt"

#ml R # commented out module load command
script=$codedir/"lift_snps.R" #restrict to eQTLs that match KGP SNPs with MAF>0.01
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile

#LD clumping==========

clumpdir=$outdir/"clump"
mkdir $clumpdir

cd $clumpdir
rm *

infile=$outdir/$tissue".signif_variant_gene_pairs.hg19.txt"
outfile=$clumpdir/$tissue".signif_variant_gene_pairs.hg19.txt"
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

ldfile=$lddir/"1000G.EUR.QC."$chr
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


