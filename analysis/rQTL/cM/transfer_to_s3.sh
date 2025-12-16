#!/bin/bash

aws s3 cp \
    dat/processed/cm_il1b/reqtl/atac.h5ad \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    dat/processed/cm_il1b/reqtl/rna.h5ad \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/sample_ids.txt \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/gene_tss_loc.bed \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/peak_loc.bed \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/pc1.txt \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 sync \
    ./dat/processed/cm_il1b/reqtl/vcfs/ \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/vcfs/

aws s3 sync \
    ./dat/processed/cm_il1b/reqtl/eqtls/ \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/eqtls/

aws s3 sync \
    ./dat/processed/cm_il1b/reqtl/caqtls/ \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/caqtls/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/cm.eqtl.pc1_reqtl.txt \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/

aws s3 cp \
    ./dat/processed/cm_il1b/reqtl/cm.eqtl.il1b_reqtl.txt \
    s3://gracie-icc/production_data/analysis_outputs/FINAL/rQTL/cM/
