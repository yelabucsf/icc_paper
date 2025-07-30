import genetics_plotting_utils as gu
import gc
import optparse

#/data/processed/QTL_calling/atac/permissive/L1_manual/B/

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

L1_colors = {
   "B": "#2D6A3F",
   "Plasma": "#7FB685",
   "CD4T": "#6B1C1C",
   "CD8T": "#A13030",
   "gdT": "#E3735E",
   "MAIT": "#FFA07A",
   "NK": "#F47024",
   "cDC": "#284F8C",
   "pDC": "#4A75D3",
   "cM": "#3E8589",
   "ncM": "#65B0B4",
   "Progen": "#595959",
   "Prolif": "#907F66"
}
L2_colors = {
   'CD4T_Naive':   '#5A1C1C',  # deep brownish-maroon
   'CD4T_Memory':  '#752211',  # deep brick red (darker, less orange)
   'Treg':         '#C45F30',  # lighter burnt orange-brown (more rusty)
   'CD8T_Naive':   '#730000',  # wine red
   'CD8T_Memory':  '#C0392B',  # classic tomato red
   'gdT':          '#E1573A',  # warm coral red
   'MAIT':         '#FF8C42',  # vivid orange
   'NK_bright':    '#F48C06',  # golden orange
   'NK_dim':       '#D95F02',  # reddish orange
   'B_Memory': '#0B5F2F', 'B_Naive': '#2D9A57', 'B_atypical': '#40BD5C', 'Plasma': '#7FC685',
   'cDC1': '#16376A', 'cDC2': '#2B5CCD', 'pDC': '#4A85E3', 'cM': '#3E8599', 'ncM': '#65B0C4', 'Progen': '#595959', 'Prolif': '#9F7F66'
}

#default is black
cell_colors={'ICC': '#000000'}
if level == 'L1_manual':
    cell_colors=L1_colors
if level == 'L2_manual':
    cell_colors=L2_colors

#cells=['B','CD4T','CD8T','gdT','MAIT','NK','cDC','cM','ncM']
for cell in cells:
    print(cell)
    print('combine chrs atac')
    atac_qtls=gu.atac_combine_chrs_with_annots(cell,0.05,str(dir+'/atac/permissive/'+level+'/'+cell+'/cis_250kb/'),str(dir+'rna/'+level+'/'+cell+'/'))#,start=5,end=6)
    print('plot atac manhattan')
    gu.plot_manhattan_custom_colors(atac_qtls ,"chrom", "pos", "pval_nominal", "eqtl_snp",outname=str(level+'_atac_'+cell+'_colored_eqtl_manhattanplot_cis250k.png'),cut=1,x_size=14, y_size=3,y_max=140,cell_hex=cell_colors[cell])
    print('plot enrich qq atac')
    gu.make_enrich_qqplot(atac_qtls, 'pval_nominal', annot_col='eqtl_snp', title="QQ Plot", cell_hex=cell_colors[cell],sample_fraction=0.5,outname=str(level+'_atac_'+cell+'_colored_eqtl_enrichqq_cis250k.png'))
    
    ##reset qtls
    #print('combine chrs rna')
    #rna_qtls=gu.rna_combine_chrs_with_annots(cell,0.05,str(dir+'/rna/'+level+'/'),str(dir+'/atac/permissive/'+level+'/'))#,start=5,end=6)
    #print('plot rna manhattan')
    #gu.plot_manhattan_custom_colors(rna_qtls ,"chrom", "pos", "pval_nominal", "atac_snp",outname=str(level+'_rna_'+cell+'_colored_caqtl_manhattanplot.png'),cut=1,x_size=14, y_size=3,y_max=120,cell_hex=cell_colors[cell])
    #print('plot enrich qq rna')
    #gu.make_enrich_qqplot(rna_qtls, 'pval_nominal', annot_col='atac_snp', title="QQ Plot", cell_hex=cell_colors[cell],sample_fraction=0.5,outname=str(level+'_atac_'+cell+'_colored_eqtl_enrichqq.png'))



print('done')
