import tempfile
import numpy as np
import pandas as pd
import coolbox
from coolbox.api import *


def plot_locus(gene_name, qtls, cell_types, out_file,
               region_size=25_000, peak=None, bed_file=None, bigwig_dir=None,
               gene_locs_file="/icc_data/References/dedup_filtered_UCSC.bed",
               gtf="/icc_data/References/gencode.v43.basic.annotation.gtf.bgz"):
    
    # read in gene locs
    gene_locs = pd.read_csv(gene_locs_file, sep="\t", names=["chr", "start", "end", "gene", "."])
    chr = gene_locs[gene_locs["gene"] == gene_name].iloc[0]["chr"]
    start = gene_locs[gene_locs["gene"] == gene_name].iloc[0]["start"]
    end = gene_locs[gene_locs["gene"] == gene_name].iloc[0]["end"]
    plot_locus = f"{chr}:{start-region_size//2}-{end+region_size//2}"

    # plot
    frame = XAxis() + GTF(gtf, row_filter="feature == 'gene'", name_attr='gene_name')
    
    if bed_file is not None:
        frame += BED(bed_file, display="collapsed", interval_height=20, labels=False)

    for i, cell_type in enumerate(cell_types):
        frame += Spacer(0.2)

        # plot eQTLs
        with tempfile.NamedTemporaryFile() as tmp:
            tmp_df = qtls["eQTLs"][cell_type].loc[(qtls["eQTLs"][cell_type]["gene"] == gene_name) &
                                                  (qtls["eQTLs"][cell_type]["snppos"] > start) &
                                                  (qtls["eQTLs"][cell_type]["snppos"] < end), 
                                                  ["snpchr", "snppos", "p-value"]]
            max_value = (-np.log10(tmp_df["p-value"])).max() + 1

            if len(tmp_df) > 0:
                tmp_df.to_csv(tmp.name, sep="\t", header=False, index=False)
                frame += Spacer(0.2)
                frame += SNP(tmp.name, col_chrom=0, col_pos=1, col_pval=2, height=2, max_value=max_value) + Title(f"{gene_name}\n{cell_type} eQTL")

        # plot caQTLs
        if peak is not None:
            with tempfile.NamedTemporaryFile() as tmp:   
                tmp_df = qtls["caQTLs"][cell_type].loc[(qtls["caQTLs"][cell_type]["gene"] == peak)  &
                                                    (qtls["caQTLs"][cell_type]["snppos"] > locus_start) &
                                                    (qtls["caQTLs"][cell_type]["snppos"] < locus_end), 
                                                    ["snpchr", "snppos", "p-value"]]
                max_value = (-np.log10(tmp_df["p-value"])).max() + 1

                if len(tmp_df) > 0:
                    tmp_df.to_csv(tmp.name, sep="\t", header=False, index=False)
                    frame += Spacer(0.2)
                    frame += SNP(tmp.name, col_chrom=0, col_pos=1, col_pval=2, height=2, max_value=max_value, threshold=0, threshold_color="purple") + Title(f"{peak}\n{cell_type} caQTL")

        # plot ATAC bigwig
        if bigwig_dir is not None:
            frame += Spacer(0.2)
            frame += BigWig(f"{bigwig_dir}/{cell_type}-TileSize-100-normMethod-ReadsInTSS.bw") + Color("black") + Title(cell_type)

    fig = frame.plot(plot_locus) 
    fig.savefig(f"{out_file}.pdf", dpi=300)
    fig.savefig(f"{out_file}.png", dpi=300)
