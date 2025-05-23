import gwaslab as gl
import optparse

#cells=['B','CD4T','CD8T','gdT','MAIT','NK','cDC','cM','ncM']
#atac_dir="/gstore/data/ctgbioinfo/gordom14/ICC/data/QTL_calling/atac/permissive/L1_manual/"
#rna_dir="/gstore/data/ctgbioinfo/gordom14/ICC/data/QTL_calling/rna/L1_manual/"

### read in data
parser = optparse.OptionParser()

parser.add_option("-c", "--cell",
                  action = "store",
                  dest = "cell", default = True,
                  help = "plink genotype file")
parser.add_option("-l", "--level",
                  action = "store",
                  dest = "level", default = True,
                  help = "prefix name")
parser.add_option("-d", "--dir",
                  action = "store",
                  dest = "dir", default = '/data2/processed/QTL_calling/',
                  help = "directory")
(options, args) = parser.parse_args()
#dir='/data2/processed/QTL_calling/'

print(options.cell)
print(options.level)
#print(options.palette)
print(options.dir)
cells=[options.cell]
level=options.level
dir=options.dir

print(cells[0])
print(level)

for cell in cells:
    print(cell)
    print(str(dir+'/atac/permissive/'+level+'/'+cell+'/cis_250kb/'+cell+'.cis_qtl_snp.@.annot.txt'))
    atac_sumstats=gl.Sumstats(str(dir+'/atac/permissive/'+level+'/'+cell+'/cis_250kb/'+cell+'.cis_qtl_snp.@.annot.txt'),fmt=None,
            rsid='variant_id',
            chrom='chrom',
            pos='pos',
            #ea='alt',
            #nea='ref',
            beta='slope',
            p='pval_nominal',
            #af='af',
            se='slope_se') 
    print('plot atac')
    atac_sumstats.plot_mqq(skip=2,save=str(cell+'_atac_mqq.png'),title=cell,)

    atac_sumstats=None

    print(str(dir+'/rna/'+level+'/'+cell+'/'+cell+'.cis_qtl_snp.@.annot.txt'))
    rna_sumstats=gl.Sumstats(str(dir+'/rna/'+level+'/'+cell+'/'+cell+'.cis_qtl_snp.@.annot.txt'),fmt=None,
            rsid='variant_id',
            chrom='chrom',
            pos='pos',
            #ea='alt',
            #nea='ref',
            beta='slope',
            p='pval_nominal',
            #af='af',
            se='slope_se') 
    print('plot rna')
    
    rna_sumstats.plot_mqq(skip=2,save=str(cell+'_rna_mqq.png'),title=cell)

print('done')
