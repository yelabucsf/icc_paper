import pandas as pd

import genetics_plotting_utils as gu

def get_AF_from_plink(freq_file,col_name='AF'):
    """
    Get the allele frequencies from the plink frequency file.
    freq_file: str
        Path to the plink frequency file.
    col_name: str
        Name of the column containing the allele frequencies.
    returns: pd.DataFrame with AF information that can be merged with the tensorQTL output.
    """
    # Load the frequency data
    freq_data = pd.read_csv(freq_file, sep='\t')

    # Extract the relevant columns
    freq_data = freq_data[['#CHROM','ID','REF','ALT','ALT_FREQS']]

    # Rename the columns for clarity
    freq_data.columns = ['chrom','variant_id','ref','alt', col_name]

    return freq_data

def make_AF_merged_file(AF_file_list, AF_names, outfile):
    # Load the raw tensorQTL data in chunks with chunksize parameter
    if AF_file_list!=None:
        # Load the allele frequency data
        i2=0
        for AF_file in AF_file_list:
            # Get the allele frequency data from the plink file
            print(AF_file)
            if i2==0:
                AF_data = get_AF_from_plink(AF_file,AF_names[i2])
                AF_data=AF_data[AF_data['variant_id']!='.']
                #AF_data=AF_data[AF_data[AF_names[i2]]!=0]
            else:
                AF_TMP=get_AF_from_plink(AF_file,AF_names[i2])
                AF_TMP=AF_TMP[AF_TMP['variant_id']!='.']
                #AF_TMP=AF_TMP[AF_TMP[AF_names[i2]]!=0]
                AF_data = pd.merge(AF_data,AF_TMP, on=['chrom','variant_id','ref','alt'], how='outer')
            i2+=1
    # filter out any sites were there are no ALT alleles for any of the data
    AF_data=AF_data[AF_data[AF_names].sum(axis=1) !=0]
    AF_data.to_csv(outfile, sep='\t', index=False) 

for i in range(1, 23):
    make_AF_merged_file([str('/data/processed/genotypes/plink/chr'+str(i)+'_ICC_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_AFR_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_AMR_freq.afreq'), str('/data/processed/genotypes/plink/chr'+str(i)+'_EAS_freq.afreq'),str('/data/processed/genotypes/plink/chr'+str(i)+'_EUR_freq.afreq')], ['ICC_AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF'],str('/data/processed/genotypes/plink/chr'+str(i)+'_merged_freqs.afreq'))



