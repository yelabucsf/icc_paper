
#import genetics_plotting_utils as gu
import gc
import optparse
import pandas as pd

#/data/processed/QTL_calling/atac/permissive/L1_manual/B/

### read in data
parser = optparse.OptionParser()

parser.add_option("-f", "--file",
                  action = "store",
                  dest = "qtl_file", default = True,
                  help = "plink genotype file")
parser.add_option("-c", "--cis",
                  action = "store",
                  dest = "cis", default = 0,
                  help = "cis window in bp default is 0")
parser.add_option("-o", "--outfile",
                  action = "store",
                  dest = "outfile", default = 'outfile',
                  help = "output file")

(options, args) = parser.parse_args()

#options needed
#qtl file
#cis
#output file
print(options.qtl_file)
print(options.cis)
print(options.outfile)

atac_qtls = pd.read_csv(options.qtl_file, sep='\t')

#add in peak location
atac_qtls['peak_start']=[int(x.replace('-',':').split(':')[1]) for x in atac_qtls['phenotype_id']]
atac_qtls['peak_end']=[int(x.replace('-',':').split(':')[2]) for x in atac_qtls['phenotype_id']]

#filter to only peaks within cis region of peak
atac_qtls_cis=atac_qtls[ (atac_qtls['peak_start'] > atac_qtls['pos']-int(options.cis) ) & (atac_qtls['peak_end'] < atac_qtls['pos']+int(options.cis) )]

atac_qtls_cis.to_csv(options.outfile, sep='\t')

##filter to only significant hits
#atac_qtls_sig=atac_qtls[ (atac_qtls['pval_nominal'] <= atac_qtls['pval_nominal_threshold']) ]
#         
#atac_qtls_sig.to_csv(options.outfile_sig, sep='\t')
#
