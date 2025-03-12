import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post
import rfunc
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")
import optparse

def is_float(element: any) -> bool:
    #If you expect None to be passed:
    if element is None: 
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False




#######################
### MAIN 
### Author Gracie Gordon 2025

##CMD example: python run_tensorqtl.py -g test_data/GEUVADIS_testdata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18 -p test_data/GEUVADIS_testdata/GEUVADIS.445_samples.expression.bed.gz -c test_data/GEUVADIS_testdata/GEUVADIS.445_samples.covariates.txt -n GEUVADIS.445_samples



### read in data
parser = optparse.OptionParser()

parser.add_option("-g", "--geno",
                  action = "store", 
                  dest = "geno", default = True,
                  help = "plink genotype file")
parser.add_option("-p", "--pheno",
                  action = "store", 
                  dest = "pheno", default = True,
                  help = "expression file with column names, chr,start,stop,gene,donor_1, .., donor_n") 
parser.add_option("-c", "--cov",
                  action = "store", 
                  dest = "cov", default = True,
                  help = "expression file with column names, chr,start,stop,gene,donor_1, .., donor_n") 
parser.add_option("-n", "--name",
                  action = "store", 
                  dest = "name", default = True,
                  help = "prefix name") 
parser.add_option("-s", "--chrom",
                  action = "store", 
                  dest = "chrom", default = 'all',
                  help = "specify chromosome to run like chr18, default is to run all")  
parser.add_option("-w", "--window",
                  action = "store", 
                  dest = "window", default = 1000000,
                  help = "cis window +/- bp, default is 1MB") 
parser.add_option("-m", "--maf",
                  action = "store", 
                  dest = "maf", default = 0.1,
                  help = "min minor allele frequency, default is 0.1") 
(options, args) = parser.parse_args()


print(str('genotypes: '+options.geno))
print(str('phenotypes: '+options.pheno))
print(str('covariates: '+options.cov))
print(str('prefix: '+options.name))
print(str('cis window: '+str(options.window)))
print(str('minor allele freq threshold: '+str(options.maf)))
print(str('chromosome: '+str(options.chrom)))

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(options.pheno)
covariates_df = pd.read_csv(options.cov, sep='\t', index_col=0).T
# convert strings to factors
#for col in covariates_df.select_dtypes(include='object').columns:

for col in covariates_df.columns:
    if is_float(covariates_df[col].iloc[0]): 
        covariates_df[col]=covariates_df[col].astype(float)
    else:
        covariates_df[col] = covariates_df[col].astype('category').cat.codes

print(covariates_df)

# PLINK reader for genotypes
pgr = pgen.PgenReader(options.geno)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

print(genotype_df.head())

# map all cis-associations (results for each chromosome are written to file)

# all genes
# cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df)

if(options.chrom == 'all'):
    # run all chrs
    print('running all chromosomes')
    cis.map_nominal(genotype_df, variant_df,
                    phenotype_df,
                    phenotype_pos_df,
                    options.name, covariates_df=covariates_df,window=int(options.window), maf_threshold=float(options.maf))


    pairs_df = pd.read_parquet(f'{options.name}.cis_qtl_pairs.{options.chrom}.parquet')

    # all genes
    # cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df)

    # genes on chr18
    cis_df = cis.map_cis(genotype_df, variant_df, 
                        phenotype_df,
                        phenotype_pos_df,
                        covariates_df=covariates_df,window=int(options.window), maf_threshold=float(options.maf), seed=123456)


    post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)


    print(cis_df.head())

    # write outputs
    print('save outputs')
    pairs_df.to_csv(f'{options.name}.cis_qtl_snp.{options.chrom}.txt', sep='\t', index=False)
    cis_df.to_csv(f'{options.name}.cis_qtl_egene.{options.chrom}.txt', sep='\t', index=False) 
else:
    print('running chromosome '+str(options.chrom))
    print(str('gtype dims: '+str(genotype_df.shape)))
    print(str('variant dims: '+str(variant_df.shape)))
    print(str('phen_df dims: '+str(phenotype_df.loc[phenotype_pos_df['chr'] == options.chrom].shape)))
    print(str('phen_df pos dims: '+str(phenotype_pos_df.loc[phenotype_pos_df['chr'] == options.chrom].shape)))
    print(str('cov dims: '+str(covariates_df.shape)))
    print('run nominal')
    cis.map_nominal(genotype_df, variant_df,
                    phenotype_df.loc[phenotype_pos_df['chr'] == options.chrom],
                    phenotype_pos_df.loc[phenotype_pos_df['chr'] == options.chrom],
                    options.name, covariates_df=covariates_df,window=int(options.window), maf_threshold=float(options.maf))


    pairs_df = pd.read_parquet(f'{options.name}.cis_qtl_pairs.{options.chrom}.parquet')

    # all genes
    # cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df)
    print('map cis')
    cis_df = cis.map_cis(genotype_df, variant_df, 
                        phenotype_df.loc[phenotype_pos_df['chr'] == options.chrom],
                        phenotype_pos_df.loc[phenotype_pos_df['chr'] == options.chrom],
                        covariates_df=covariates_df,window=int(options.window), maf_threshold=float(options.maf), seed=123456)
    print('calc qvals')

    post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)


    print(cis_df.head())

    # write outputs
    print('save outputs')
    pairs_df.to_csv(f'{options.name}.cis_qtl_snp.{options.chrom}.txt', sep='\t', index=False)
    cis_df.to_csv(f'{options.name}.cis_qtl_egene.{options.chrom}.txt', sep='\t', index=True) 

