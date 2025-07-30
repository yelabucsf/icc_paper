
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
                  help = "qtl snp file")
parser.add_option("-e", "--egenefile",
                  action = "store",
                  dest = "egene_file", default = True,
                  help = "qtl egene file")
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

qtls = pd.read_csv(options.qtl_file, sep='\t')
egenes = pd.read_csv(options.egene_file, sep='\t')

#merge
qtls = qtls.merge(egenes[['phenotype_id','pvalue_nominal_threshold']], how='left', on='phenotype_id')
##filter to only significant hits
qtls_sig=qtls[ (qtls['pval_nominal'] <= qtls['pval_nominal_threshold']) ]
      
#save
qtls_sig.to_csv(options.outfile, sep='\t')

