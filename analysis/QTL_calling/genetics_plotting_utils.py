import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import scanpy as sc
from tensorqtl import pgen
from scipy.stats import probplot
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare
import gwaslab as gl
from scipy.stats import probplot
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare
from tqdm import tqdm
import dask
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

def get_gene_id_mapping_dask(adata_file):

    """
    Get the gene id mapping from the AnnData object.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    """
    # Load the AnnData object
    adata = sc.read_h5ad(adata_file)

    # Create a DataFrame with gene_id and phenotype_id
    gene_id_mapping = pd.DataFrame(adata.var['name'])
    gene_id_mapping['phenotype_id'] = gene_id_mapping.index
    gene_id_mapping.reset_index(drop=True, inplace=True)
    gene_id_mapping.columns = ['gene_id', 'phenotype_id']
    adata=0
    return dd.from_pandas(gene_id_mapping)

def get_AF_from_plink_dask(freq_file,col_name='AF'):
    """
    Get the allele frequencies from the plink frequency file.
    freq_file: str
        Path to the plink frequency file.
    col_name: str
        Name of the column containing the allele frequencies.
    returns: pd.DataFrame with AF information that can be merged with the tensorQTL output.
    """
    # Load the frequency data
    freq_data = dd.read_csv(freq_file, sep='\t',dtype={'ALT_FREQS': 'float64'})
    
    # Extract the relevant columns
    freq_data = freq_data[['#CHROM','ID','REF','ALT','ALT_FREQS']]
    
    # Rename the columns for clarity
    freq_data.columns = ['chrom','variant_id','ref','alt', col_name]
    
    return freq_data


def save_dask_dataframe_to_single_file(dask_dataframe, output_file):
    """
    Save a Dask DataFrame to a single file by appending each partition.

    Parameters:
    - dask_dataframe: Dask DataFrame to save.
    - output_file: str, the file where the data will be saved.
    """
    # Ensure the output directory exists
    dask.config.set(temporary_directory='/data2/tmp')
    import os
    import gc
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Iterate over partitions and append each one to the file
    with ProgressBar():
        for i, partition in enumerate(dask_dataframe.to_delayed()):
            mode = 'w' if i == 0 else 'a'  # Write mode for the first partition, append mode for others
            header = (i == 0)  # Write header only for the first partition
            print(i)
            partition.compute().to_csv(output_file, mode=mode, index=False, header=header)
            del partition
            gc.collect()

def process_raw_tensorqtl_dask(adata_file, geno_pvar_file, tensor_out_file, outfile_name="None", AF_file_list=None, AF_names=None,block_size="1MB"):
    """
    Process the raw tensorQTL data using Dask for improved performance.

    Parameters:
    - adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format.
    - geno_pvar_file: str
        Path to the .pvar file containing the genotype data.
    - tensor_out_file: str
        Path to the tensorQTL output file.
    - outfile_name: str
        Name of the output file to save the processed data.
    - AF_file_list: list
        List of allele frequency files.
    - AF_names: list
        List of column names for allele frequencies.

    Returns:
    - If outfile_name is not provided, returns a Dask DataFrame with the processed data.
    """
    # Load SNP data
    print('load snps')
    genos_pos = dd.read_csv(geno_pvar_file, skiprows=14, sep='\t')
    genos_pos.columns = ['chrom', 'pos', 'variant_id', 'ref', 'alt', 'info']

    # Load gene ID info
    print('load genes')
    gene_id_mapping = get_gene_id_mapping_dask(adata_file)

    print('load AF')
    # Load allele frequency data if provided
    if AF_file_list is not None:
        i2 = 0
        for AF_file in AF_file_list:
            if i2 == 0:
                AF_data = get_AF_from_plink_dask(AF_file, AF_names[i2])
                AF_data = AF_data[AF_data['variant_id'] != '.']
            else:
                AF_TMP = get_AF_from_plink_dask(AF_file, AF_names[i2])
                AF_TMP = AF_TMP[AF_TMP['variant_id'] != '.']
                AF_data = dd.merge(AF_data, AF_TMP, on=['chrom', 'variant_id', 'ref', 'alt'], how='outer')
            i2 += 1

    print('load tensorQTL')
    # Load the tensorQTL data using Dask
    tensor_data = dd.read_csv(tensor_out_file, sep='\t', blocksize=block_size)
    print(tensor_data.npartitions)
    # Annotate the data
    print('merge gene_id')
    tensor_data = dd.merge(tensor_data,gene_id_mapping, on='phenotype_id', how='left')
    print('merge snp_id')
    tensor_data = dd.merge(tensor_data,genos_pos, on='variant_id', how='left')

    print('merge AF')
    if AF_file_list is not None:
        tensor_data = tensor_data.merge(AF_data, on=['chrom', 'variant_id', 'ref', 'alt'], how='left')
    #print('repartition')
    # repartition
    #tensor_data=tensor_data.repartition(partition_size=block_size)
    print(tensor_data.npartitions)
    # Save the annotated data to a file or return it
    print('save files')
    if outfile_name != "None":
        # Save the Dask DataFrame in chunks
        save_dask_dataframe_to_single_file(tensor_data, outfile_name)
        #with ProgressBar():
        #    tensor_data.to_csv(outfile_name, sep='\t', index=False, single_file=True)
    else:
        return tensor_data

def get_gene_id_mapping(adata_file):

    """
    Get the gene id mapping from the AnnData object.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    """
    # Load the AnnData object
    adata = sc.read_h5ad(adata_file)

    # Create a DataFrame with gene_id and phenotype_id
    gene_id_mapping = pd.DataFrame(adata.var['name'])
    gene_id_mapping['phenotype_id'] = gene_id_mapping.index
    gene_id_mapping.reset_index(drop=True, inplace=True)
    gene_id_mapping.columns = ['gene_id', 'phenotype_id']
    adata=0
    return gene_id_mapping

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

def process_raw_tensorqtl_new(gene_id_file,geno_pvar_file, tensor_out_file, outfile_name="None", chunk_size='None', AF_file='None'):
    """
    Process the raw tensorQTL data and return a DataFrame with the relevant information added.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    geno_pvar_file: str
        Path to the .pvar file containing the genotype data, which includes the variant IDs and their positions.
    tensor_out_file: str
        Path to the tensorQTL output file, which contains the raw QTL results.

    Note: If chunk_size is used, the file will be saved to a file and not returned
    """
    import os
    ## skip if file does not exist
    if not os.path.exists(tensor_out_file):
        print('file does not exist')
        return

    # Method 1: Check if the directory exists first
    path = os.path.dirname(outfile_name)
    print(path)
    if not os.path.exists(path):
        os.mkdir(path)

    # load gene id info
    #gene_id_mapping=get_gene_id_mapping(adata_file)
    gene_id_mapping=pd.read_csv(gene_id_file,sep='\t')
    #print(gene_id_mapping.head())
    if AF_file != 'None':
        AF_data=pd.read_csv(AF_file,sep='\t')
        #print(AF_data.head())

    # Load snp data
    genos_pos=pd.read_csv(geno_pvar_file,skiprows=14,sep='\t')
    genos_pos.columns=['chrom','pos','variant_id','ref','alt','info']
    #print(genos_pos.shape)
    genos_pos=genos_pos[genos_pos['variant_id'].isin(AF_data['variant_id'])]
    #print(genos_pos.shape)
    
    # Load the raw tensorQTL data in chunks with chunksize parameter
    print('begin merge')
    if chunk_size!='None':
        i=1
        # Get the total number of rows in the CSV for progress tracking
        total_rows = sum(1 for _ in open(tensor_out_file)) - 1  # Subtract 1 for the header row

        # Initialize tqdm progress bar
        tqdm_bar = tqdm(total=total_rows, desc="Processing CSV")
        for chunk in pd.read_csv(tensor_out_file, sep='\t', chunksize=chunk_size):
            # annotate
            
            chunk_annot=pd.merge(genos_pos,chunk,on='variant_id', how='right')
            #chunk_annot=pd.merge(chunk_annot,gene_id_mapping,on='phenotype_id', how='left')
            chunk_annot=pd.merge(chunk_annot,gene_id_mapping,on='phenotype_id', how='left')
            print(chunk_annot)
            if AF_file!=None:
                chunk_annot=pd.merge(chunk_annot,AF_data, on=['chrom','variant_id','ref','alt'], how='left')
                #chunk_annot=pd.merge(AF_data, chunk_annot, on=['variant_id'], how='left')
             
            chunk_annot=chunk_annot[chunk_annot['variant_id'] !='.']
            print(chunk_annot.loc[chunk_annot['chrom'].isna()])


            #fix integers
            chunk_annot['chrom']=chunk_annot['chrom'].astype(int)
            chunk_annot['pos']=chunk_annot['pos'].astype(int)
            if(i==1):
                chunk_annot.to_csv(outfile_name, sep='\t', index=False)
                i+=1

            else:
                chunk_annot.to_csv(outfile_name, sep='\t', index=False, header=False, mode='a')
            tqdm_bar.update(len(chunk))
    else:
        # run without streaming data for small files
        qtl_snp_output=pd.read_csv(tensor_out_file, sep='\t')

        # annotate
        #qtl_snp_output_annot=pd.merge(qtl_snp_output,genos_pos,on='variant_id', how='left')
        qtl_snp_output_annot=pd.merge(qtl_snp_output_annot,gene_id_mapping,on='phenotype_id', how='left')

        if AF_file!=None:
            qtl_snp_output_annot=pd.merge(qtl_snp_output_annot,AF_data, on=['chrom','variant_id','ref','alt'], how='left')




        if(outfile_name!="None"):
            qtl_snp_output_annot.to_csv(outfile_name, sep='\t', index=False)

        # Return the processed data
        return qtl_snp_output_annot


def process_raw_tensorqtl(adata_file, geno_pvar_file,tensor_out_file, outfile_name="None", chunk_size='None', AF_file_list=None,AF_names=None):
    """
    Process the raw tensorQTL data and return a DataFrame with the relevant information added.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    geno_pvar_file: str
        Path to the .pvar file containing the genotype data, which includes the variant IDs and their positions.
    tensor_out_file: str
        Path to the tensorQTL output file, which contains the raw QTL results.
    
    Note: If chunk_size is used, the file will be saved to a file and not returned
    """
    ## skip if file does not exist
    if not os.path.exists(tensor_out_file):
        print('file does not exist')
        return

    # Load snp data
    genos_pos=pd.read_csv(geno_pvar_file,skiprows=14,sep='\t')
    genos_pos.columns=['chrom','pos','variant_id','ref','alt','info']

    # load gene id info
    gene_id_mapping=get_gene_id_mapping(adata_file)

    # Load the raw tensorQTL data in chunks with chunksize parameter
    if AF_file_list!=None:
        # Load the allele frequency data
        i2=0
        for AF_file in AF_file_list:
            # Get the allele frequency data from the plink file
            #print(AF_file)
            if i2==0:
                AF_data = get_AF_from_plink(AF_file,AF_names[i2])
                AF_data=AF_data[AF_data['variant_id']!='.']
            else:
                AF_TMP=get_AF_from_plink(AF_file,AF_names[i2])
                AF_TMP=AF_TMP[AF_TMP['variant_id']!='.']
                AF_data = pd.merge(AF_data,AF_TMP, on=['chrom','variant_id','ref','alt'], how='outer')
            i2+=1
  
    if chunk_size!='None':
        i=1
        # Get the total number of rows in the CSV for progress tracking
        total_rows = sum(1 for _ in open(tensor_out_file)) - 1  # Subtract 1 for the header row

        # Initialize tqdm progress bar
        tqdm_bar = tqdm(total=total_rows, desc="Processing CSV")
        for chunk in pd.read_csv(tensor_out_file, sep='\t', chunksize=chunk_size):
            # annotate
            chunk_annot=pd.merge(chunk,genos_pos,on='variant_id', how='left')
            chunk_annot=pd.merge(chunk_annot,gene_id_mapping,on='phenotype_id', how='left')

            if AF_file_list!=None:
                chunk_annot=pd.merge(chunk_annot,AF_data, on=['chrom','variant_id','ref','alt'], how='left')
                
            if(i==1):
                chunk_annot.to_csv(outfile_name, sep='\t', index=False)
                i+=1
                
            else:
                chunk_annot.to_csv(outfile_name, sep='\t', index=False, header=False, mode='a')
            tqdm_bar.update(len(chunk))

    else:
        # run without streaming data for small files
        qtl_snp_output=pd.read_csv(tensor_out_file, sep='\t')

        # annotate
        qtl_snp_output_annot=pd.merge(qtl_snp_output,genos_pos,on='variant_id', how='left')
        qtl_snp_output_annot=pd.merge(qtl_snp_output_annot,gene_id_mapping,on='phenotype_id', how='left')
        
        if AF_file_list!=None:
            qtl_snp_output_annot=pd.merge(qtl_snp_output_annot,AF_data, on=['chrom','variant_id','ref','alt'], how='left')




        if(outfile_name!="None"):
            qtl_snp_output_annot.to_csv(outfile_name, sep='\t', index=False)
        
        # Return the processed data
        return qtl_snp_output_annot


def process_raw_tensorqtl_old(adata_file, geno_pvar_file,tensor_out_file, outfile_name="None"):
    """
    Process the raw tensorQTL data and return a DataFrame with the relevant information added.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    geno_pvar_file: str
        Path to the .pvar file containing the genotype data, which includes the variant IDs and their positions.
    tensor_out_file: str
        Path to the tensorQTL output file, which contains the raw QTL results.
    """
    # Load the raw tensorQTL data
    qtl_snp_output=pd.read_csv(tensor_out_file, sep='\t')

    # Load snp data
    genos_pos=pd.read_csv(geno_pvar_file,skiprows=14,sep='\t')
    genos_pos.columns=['chrom','pos','variant_id','ref','alt','info']

    # load gene id info
    gene_id_mapping=get_gene_id_mapping(adata_file)

    # annotate
    qtl_snp_output_annot=pd.merge(qtl_snp_output,genos_pos,on='variant_id', how='left')
    qtl_snp_output_annot=pd.merge(qtl_snp_output_annot,gene_id_mapping,on='phenotype_id', how='left')

    if(outfile_name!="None"):
        qtl_snp_output_annot.to_csv(outfile_name, sep='\t', index=False)
    
    # Return the processed data
    return qtl_snp_output_annot


def merge_per_chr_annotated_tensorqtl(adata_file,genotype_dir,qtl_dir,celltype):
    """
    Merge the per chromosome annotated tensorQTL results into a single DataFrame.
    
    """
    # Initialize an empty DataFrame to store the merged results
    merged_df = pd.DataFrame()

    # Loop through each chromosome and read the corresponding file
    for i in range(1, 23):
        #file_path = f"chr{i}_annotated_tensorqtl.txt"
        df = process_raw_tensorqtl(adata_file,
                               str(genotype_dir+'chr'+str(i)+'_cleaned.pvar'),
                               str(qtl_dir+celltype+'/'+celltype+'.cis_qtl_snp.'+str(i)+'.txt'))
        merged_df = pd.concat([merged_df, df], ignore_index=True)
    return(merged_df)


def calc_lambdagc(data_original, column):
    # Calcluate lambda GC on full data
    chi2_statistic, p_value = chisquare(data_original[column])
    denom,pvalue=chisquare(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_original[column]))))
    lambdagc = np.median(chi2_statistic)/denom
    return(lambdagc)


from scipy.stats import probplot
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare

def calc_lambdagc(data_original, column):
    # Calcluate lambda GC on full data
    chi2_statistic, p_value = chisquare(data_original[column])
    denom,pvalue=chisquare(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_original[column]))))
    lambdagc = np.median(chi2_statistic)/denom
    return(lambdagc)


def make_qqplot(data_original, column, title="QQ Plot", sample_fraction=1):
    """
    Generate a QQ plot for a specified column in a DataFrame and calculate lambdaGC.

    Parameters:
    - data: pandas DataFrame containing the data.
    - column: str, the column name for which to generate the QQ plot.
    - title: str, the title of the plot (default is "QQ Plot").
    - sample_fraction: float, random fraction of data to plot (default is 1, meaning no sampling).
    """
    # Calcluate lambda GC on full data
    #chi2_statistic, p_value = chisquare(data_original[column])
    #denom,pvalue=chisquare(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_original[column]))))
    #lambdagc = np.median(chi2_statistic)/denom
    lambdagc = calc_lambdagc(data_original, column)

    data=data_original#.sample(frac=sample_fraction, random_state=1)

    plt.figure(figsize=(8, 8))
    #-log10
    obs=-1*np.log10(data[column])
    #probplot(data[column], dist="uniform", plot=plt)
    #probplot(data[column], dist=stats.loguniform,  plot=plt)
    #exp=-1*np.log10(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data[column]))))
    exp=-1*np.log10([ i/len(data[column]) for i in range(1,len(data[column])+1)])
    df = pd.DataFrame(list(zip(sorted(obs), sorted(exp))),
               columns =['obs', 'exp'])
    
    df=df.sample(frac=sample_fraction, random_state=1)
    sns.set_context("notebook", font_scale=2)
    sns.scatterplot(df, x='exp', y='obs',edgecolor=None)
    plt.axline((0, 0), slope=1, color='red', linestyle='--')


    plt.title(title+f" (lambda GC: {lambdagc:.2f})")
    plt.show()

def make_enrich_qqplot(data_original, column, annot_col='eqtl_snp', title="QQ Plot", cell_hex='None',sample_fraction=1,outname="None"):
    """
    Generate a QQ plot for a specified column in a DataFrame and calculate lambdaGC.

    Parameters:
    - data: pandas DataFrame containing the data.
    - column: str, the column name for which to generate the QQ plot.
    - title: str, the title of the plot (default is "QQ Plot").
    - sample_fraction: float, random fraction of data to plot (default is 1, meaning no sampling).
    """
    # Calcluate lambda GC on full data
    #chi2_statistic, p_value = chisquare(data_original[column])
    #denom,pvalue=chisquare(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_original[column]))))
    #lambdagc = np.median(chi2_statistic)/denom
    lambdagc = calc_lambdagc(data_original, column)

    data=data_original#.sample(frac=sample_fraction, random_state=1)

    plt.figure(figsize=(8, 8))

    #split by annotation
    data_annot1=data[data[annot_col]=='Yes']
    data_annot2=data[data[annot_col]=='No']

    #-log10
    obs=-1*np.log10(data_annot1[column])
    #probplot(data[column], dist="uniform", plot=plt)
    #probplot(data[column], dist=stats.loguniform,  plot=plt)
    #exp=-1*np.log10(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_annot1[column]))))
    exp=-1*np.log10([ i/len(data_annot1[column]) for i in range(1,len(data_annot1[column])+1)])
    
    df = pd.DataFrame(list(zip(sorted(obs), sorted(exp))),
               columns =['obs', 'exp'])
    df[annot_col]='Yes'

    obs=-1*np.log10(data_annot2[column])
    #exp=-1*np.log10(stats.uniform.ppf(np.linspace(0.001, 0.999, len(data_annot2[column]))))
    exp=-1*np.log10([ i/len(data_annot2[column]) for i in range(1,len(data_annot2[column])+1)])
    df2 = pd.DataFrame(list(zip(sorted(obs), sorted(exp))),
               columns =['obs', 'exp'])
    df2[annot_col]='No'
    df_combined = pd.concat([df, df2], ignore_index=True)
    df_combined=df_combined.sample(frac=sample_fraction, random_state=1)

    hue_order=['Yes','No']
    if cell_hex=='None':
        new_palette = {'Yes': sns.color_palette()[1],'No': sns.color_palette()[0] }  # Custom colors for the groups
    else:
        new_hex= lighten_color(cell_hex, 0.5)
        print(new_hex)
        new_palette = {'Yes': cell_hex, 'No': new_hex }  # Custom colors for the groups
    
    sns.set_context("notebook", font_scale=2)
    sns.scatterplot(df_combined, 
                    x='exp', 
                    y='obs',
                    hue=annot_col,
                    hue_order=hue_order,  # Specify the order of the colors
                    palette=new_palette,
                    edgecolor=None)
    plt.axline((0, 0), slope=1, color='red', linestyle='--')


    plt.title(title)
    if outname !="None":
        
        plt.savefig(outname, format='png',  dpi=300)
    #plt.show()



def plot_boxplot_with_points(data, x_col, y_col, title="Boxplot with Points"):
    """
    Create a boxplot with individual data points overlaid.

    Parameters:
    - data: pandas DataFrame containing the data.
    - x_col: str, the column name for the x-axis (categorical variable).
    - y_col: str, the column name for the y-axis (numerical variable).
    - title: str, the title of the plot.
    """
    plt.figure(figsize=(10, 6))
    sns.set_context("notebook", font_scale=1.5)
    
    # Create the boxplot
    sns.boxplot(x=x_col, y=y_col, data=data, showfliers=False)#, palette="Set2")
    
    # Overlay individual data points
    sns.stripplot(x=x_col, y=y_col, data=data, color="black", alpha=0.6, jitter=True)
    
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.show()



def prep_boxplot_df(qtl_output, plink_prefix, pheno_id, variant_id,adata_file, pheno_col='hgnc_id'):
    """
    Create a DataFrame for boxplotting with the relevant information added.
    qtl_output: str
        Path to the .txt file containing the QTL results.
    plink_prefix: str
        Path to the prefix of the PLINK files (without file extension).
    pheno_id: str
        The phenotype ID to filter the data.
    variant_id: str
        The variant ID to filter the data.
    adata_file: str
        Path to the .h5ad file containing the gene expression data in AnnData format, where ensmbl IDs are in index and gene name
    pheno_col: str
        The column name in the QTL output file that contains the phenotype IDs (default is 'hgnc_id').
    TODO: add ability to show Alleles tested outside of just dosage, condider adding ability to regress out covariates from phenotype file and plot the residuals. 
    """

    # Get mapping
    gene_mapping=get_gene_id_mapping(adata_file)
    gene_mapping.columns=['hgnc_id','gene_id']

    # load in phenotype info
    pheno=pd.read_csv(qtl_output, sep='\t')
    pheno_annot=pd.merge(gene_mapping,pheno,on='gene_id', how='right')

    #filter one gene
    pheno_annot_filtered=pheno_annot[pheno_annot[pheno_col]==pheno_id]
    df_pheno=pheno_annot_filtered.iloc[:,5:].T
    df_pheno.columns=[pheno_id]
    df_pheno


    # Get the SNP data
    snps=pgen.PgenReader(plink_prefix)
    cur_snp=snps.read(variant_id)

    df_pheno=pd.merge(df_pheno,cur_snp, how='inner', left_index=True, right_index=True)


    return(df_pheno)


def make_maimi_plots(snp_output_1, snp_output_2, ):
    # this takes forever >1 hour
    sumstats1=gl.Sumstats(snp_output_1,
            fmt=None,
            rsid='variant_id',
            chrom='chrom',
            pos='pos',
            ea='alt',
            nea='ref',
            beta='slope',
            p='pval_nominal',
            af='af',
            se='slope_se'
            )
    sumstats2=gl.Sumstats(snp_output_2,
                fmt=None,
                rsid='variant_id',
                chrom='chrom',
                pos='pos',
                ea='alt',
                nea='ref',
                beta='slope',
                p='pval_nominal',
                af='af',
                se='slope_se'
                )
    
    gl.plot_miami(sumstats1, sumstats2)

def gwaslab_manhattan(qtl_snp_output_annot, outname="None",cut_val=0,dpi=300):
    """
    Convert tensorQTL sumstats to a format that can be used in the QTLtools pipeline.
    """
    # import tensorqtl object previously processed with function 'process_raw_tensorqtl'
    # and convert to QTLtools format
    sumstats=gl.Sumstats(qtl_snp_output_annot,
            fmt=None,
            rsid='variant_id',
            chrom='chrom',
            pos='pos',
            ea='alt',
            nea='ref',
            beta='slope',
            p='pval_nominal',
            af='af',
            se='slope_se'
            )
    
    sumstats.plot_mqq(cut=cut_val)
    if outname!="None":
        sumstats.plot_mqq(cut=cut_val,save=outname,save_args={"dpi":dpi,"facecolor":"white"})

    
#convert_tensorqtl_sumstats(qtl_snp_output_annot)    
    

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def lighten_color(hex_color, amount=0.5):
    """Lightens a given hex color.

    Args:
        hex_color: The hex color string (e.g., "#RRGGBB").
        amount: The amount to lighten the color, between 0 and 1. 
            0 returns the original color, 1 returns white.

    Returns:
        A lighter hex color string.
    """
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    r = int(r + (255 - r) * amount)
    g = int(g + (255 - g) * amount)
    b = int(b + (255 - b) * amount)

    return "#{:02x}{:02x}{:02x}".format(r, g, b)

# Example usage
#original_color = "#3498db"
#lighter_color = lighten_color(original_color, 0.8)
#print(f"Original color: {original_color}")
#print(f"Lighter color: {lighter_color}")

def plot_manhattan_custom_colors(data, chrom_col, pos_col, pval_col, color_col,cut=0,outname='None', title="Manhattan Plot", xlabel="Chromosome", ylabel="-log10(p-value)",x_size=12, y_size=6,dpi=300, y_max='None',cell_hex='None',legend='No',flip_pal='No'):
    """
    Plot a Manhattan plot with custom colors for each point.

    Parameters:
    - data: pandas DataFrame containing the data.
    - chrom_col: str, the column name for chromosome numbers.
    - pos_col: str, the column name for base pair positions.
    - pval_col: str, the column name for p-values.
    - color_col: str, the column name for custom colors.
    - title: str, the title of the plot.
    - xlabel: str, the label for the x-axis.
    - ylabel: str, the label for the y-axis.
    - x_size: int, the width of the plot.
    - y_size: int, the height of the plot.
    - dpi: int, the resolution of the plot.
    - cell_hex: str, the hex color code for the points.
    """
    # Calculate -log10(p-value)
    data['-log10(p-value)'] = -np.log10(data[pval_col])

    # Filter out points with p-values greater than the cut-off
    if cut>0:
        data = data[data['-log10(p-value)'] > cut]

    # Sort data by chromosome and position
    data = data.sort_values(by=[chrom_col, pos_col,color_col])
    #print(data)
    # Create a new column for cumulative positions
    data['cumulative_pos'] = data.groupby(chrom_col)[pos_col].cumsum()
    data['current_position']=range(len(data['pos']))

    hue_order=['Yes','No']

    data = data.sort_values(by=[chrom_col, color_col])
    #print(data)
    # Create a mapping for chromosome midpoints for x-axis labels
    chrom_midpoints = data.groupby(chrom_col)['current_position'].mean()
    #print(chrom_midpoints)
    plt.figure(figsize=(x_size, y_size))

    #create palatte from cell hex codes 
    if cell_hex=='None':
        new_palette = {'Yes': sns.color_palette()[1],'No': sns.color_palette()[0] }  # Custom colors for the groups
    else:
        new_hex= lighten_color(cell_hex, 0.5)
        if(flip_pal=='Yes'):
            new_palette = {'Yes': new_hex,'No': cell_hex }
        else:
            new_palette = {'Yes': cell_hex, 'No': new_hex }  # Custom colors for the groups
    if(y_max=='None'):
        y_max=data['-log10(p-value)'].max() + 1
    sns.scatterplot(
        data=data,
        x='current_position',
        y='-log10(p-value)',
        hue=color_col,
        s=8,
        #palette="viridis",
        edgecolor=None,
        hue_order=hue_order,  # Specify the order of the colors
        palette=new_palette,
        legend='full'
    ).set_ylim(0,y_max)
    if legend=='Yes':
        # Customize the legend
        plt.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc='upper left')

    # Set x-axis ticks to chromosome midpoints
    plt.xlim(data["current_position"].min(), data["current_position"].max())
    plt.xticks(ticks=chrom_midpoints, labels=chrom_midpoints.index, rotation=45)


    # Add labels and title
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Add a horizontal line for genome-wide significance threshold (optional)
    plt.axhline(y=-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide significance')

    plt.legend(
        #title=color_col,
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        loc='upper left'           # Align the legend to the upper left of the anchor
    )
    plt.tight_layout()
    
    if outname !="None":
        
        plt.savefig(outname, format='png',  dpi=300)
    plt.show()


def locuszoom_plot_qtl(qtl_snp_output_annot,geneid,chr, start, stop, gtf_select='refseq_hg38_gtf',ld_select='1kg_pan_hg38', cis=0):

    """
    Plot a locuszoom plot for the QTL results.
    qtl_snp_output_annot: DataFrame
        The QTL results with the relevant information added.
    geneid: str
        The gene ID to filter the data.
    chr: int
        The chromosome number.  
    start: int
        The start position of the region to plot.
    stop: int
        The stop position of the region to plot.
    gtf_select: str
        The GTF file to use for annotation (default is 'refseq_hg38_gtf').
    ld_select: str
        The LD file to use for annotation (default is '1kg_pan_hg38').
    cis: int
        The cis region to plot (default is 0).
    Returns: plot

    Recommeneded to use a single chromosome QTL output file for this function.
    """

    ## note the reference files must be downloaded first gl.download_ref("ensembl_hg38_gtf"), gl.download_ref("1kg_pan_hg38")
    # import tensorqtl object previously processed with function 'process_raw_tensorqtl'
    # and convert to QTLtools format
    sumstats=gl.Sumstats(qtl_snp_output_annot[qtl_snp_output_annot['gene_id']==geneid],
            fmt=None,
            rsid='variant_id',
            chrom='chrom',
            pos='pos',
            ea='alt',
            nea='ref',
            beta='slope',
            p='pval_nominal',
            af='af',
            se='slope_se'
            )
    print(sumstats.basic_check())
    print(sumstats.data.head())
    sumstats.plot_mqq(mode="r", region=(chr,start-cis,stop+cis),build="38",skip=2, #anno=True, #anno_source='ensembl', 
                      vcf_path=gl.get_path(ld_select),#gtf_path=gl.get_path(gtf_select),
                      taf=[4,0,0.95,1,1])#region_recombination=False,
    #sumstats.plot_mqq(mode="r",build="38",region=(chr,start-cis,stop+cis))

def annotate_qtls_with_other_qtls(qtl_result_df_res1,qtl_result_df_res2, permutation_file_res2,annot_colname='eqtl_snp'):
    """
    annotate tensorQTL results with other QTLs results
    """
    perm_results=pd.read_csv(permutation_file_res2, sep='\t')
    #filter to only results passing nominal permutation testing
    qtl_result_df_res2=qtl_result_df_res2.merge(perm_results[perm_results['pval_nominal'] < perm_results['pval_nominal_threshold']][['phenotype_id','pval_nominal_threshold']], on='phenotype_id', how='left')
    sig_annots_snps=qtl_result_df_res2[qtl_result_df_res2['pval_nominal']<qtl_result_df_res2['pval_nominal_threshold']]['variant_id'].unique().tolist()
    qtl_result_df_res1[annot_colname]='No'

    qtl_result_df_res1.loc[qtl_result_df_res1['variant_id'].isin(sig_annots_snps), annot_colname] = 'Yes'
    return(qtl_result_df_res1)

def atac_combine_chrs_with_annots(cell, pval_thresh, atac_dir, rna_dir, annot_colname='eqtl_snp', qtl_prefix='.cis_qtl_snp.',qtl_suffix='.annot.txt',start=1, end=23):
    for i in range(start, end):
        print(str('chr:'+str(i)))
        atac_qtl_file = f"{atac_dir}/{cell}.cis_qtl_snp.{i}.annot.txt"
        rna_qtl_file = f"{rna_dir}/{cell}.cis_qtl_snp.{i}.annot.txt"
        if i == start:
            atac_qtls = pd.read_csv(atac_qtl_file, sep='\t')
            rna_qtls = pd.read_csv(rna_qtl_file, sep='\t')
            atac_qtls=annotate_qtls_with_other_qtls(atac_qtls,rna_qtls, f"{rna_dir}/{cell}.cis_qtl_egene.{i}.txt", annot_colname='eqtl_snp')
            #optional filtering of pvalues to save RAM
            atac_qtls = atac_qtls[atac_qtls['pval_nominal'] < pval_thresh]
        else:
            atac_qtls_tmp= pd.read_csv(atac_qtl_file, sep='\t')
            atac_qtls_tmp = atac_qtls_tmp[atac_qtls_tmp['pval_nominal'] < pval_thresh]
            rna_qtls_tmp = pd.read_csv(rna_qtl_file, sep='\t')
            atac_qtls_tmp=annotate_qtls_with_other_qtls(atac_qtls_tmp,rna_qtls_tmp, f"{rna_dir}/{cell}.cis_qtl_egene.{i}.txt", annot_colname='eqtl_snp')

            atac_qtls = pd.concat([atac_qtls, atac_qtls_tmp], ignore_index=True)
            #optional filtering of pvalues to save RAM
            atac_qtls = atac_qtls[atac_qtls['pval_nominal'] < pval_thresh]
    return(atac_qtls)

def rna_combine_chrs_with_annots(cell, pval_thresh, rna_dir, atac_dir, annot_colname='atac_snp', qtl_prefix='.cis_qtl_snp.',qtl_suffix='.annot.txt',start=1, end=23):
    """
    Combine chromosomes into a single dataframe for plotting.
    """
    for i in range(start, end):
        print(str('chr:'+str(i)))
        atac_qtl_file = f"{atac_dir}/{cell}.cis_qtl_snp.{i}.annot.txt"
        rna_qtl_file = f"{rna_dir}/{cell}.cis_qtl_snp.{i}.annot.txt"
        if i == start:
            atac_qtls = pd.read_csv(atac_qtl_file, sep='\t')
            atac_qtls = atac_qtls[atac_qtls['pval_nominal'] < pval_thresh]
            rna_qtls = pd.read_csv(rna_qtl_file, sep='\t')
            rna_qtls=annotate_qtls_with_other_qtls(rna_qtls,atac_qtls, f"{atac_dir}/{cell}.cis_qtl_egene.{i}.txt", annot_colname=annot_colname)
            #optional filtering of pvalues to save RAM
            rna_qtls = rna_qtls[rna_qtls['pval_nominal'] < pval_thresh]
        else:
            atac_qtls_tmp= pd.read_csv(atac_qtl_file, sep='\t')
            atac_qtls_tmp = atac_qtls_tmp[atac_qtls_tmp['pval_nominal'] < pval_thresh]
            rna_qtls_tmp = pd.read_csv(rna_qtl_file, sep='\t')
            rna_qtls_tmp=annotate_qtls_with_other_qtls(rna_qtls,atac_qtls_tmp, f"{atac_dir}/{cell}.cis_qtl_egene.{i}.txt", annot_colname='atac_snp')

            rna_qtls = pd.concat([rna_qtls, rna_qtls_tmp], ignore_index=True)
            #optional filtering of pvalues to save RAM
            rna_qtls = rna_qtls[rna_qtls['pval_nominal'] < pval_thresh]
    return(rna_qtls)










































































































































































