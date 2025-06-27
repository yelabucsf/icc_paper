#!/bin/bash

# L1 permissive
Rscript annotate_peaks.R --peak_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L1_peaks/XHLT1-POOL-XX1-SC1/peaks.bed --out_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L1_peaks/XHLT1-POOL-XX1-SC1/peaks_annotated.bed

# L1 stringent
Rscript annotate_peaks.R --peak_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L1_peaks/PeakCalls/stringent_peaks/stringent_L1_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping.bed --out_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L1_peaks/PeakCalls/stringent_peaks/stringent_L1_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping_annotated.bed

# L2_permissive
Rscript annotate_peaks.R --peak_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/XHLT1-POOL-XX1-SC1/peaks.bed --out_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/XHLT1-POOL-XX1-SC1/peaks_annotated.bed

# L2_stringent
Rscript annotate_peaks.R --peak_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/PeakCalls/stringent_peaks/stringent_L2_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping.bed --out_bed /data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/PeakCalls/stringent_peaks/stringent_L2_peaks_spm5e250_rulen12/All_Samples.fwp.filter.non_overlapping_annotated.bed