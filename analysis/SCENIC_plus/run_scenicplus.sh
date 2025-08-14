#!/bin/bash

mamba activate scenicplus
python run_pycistopic.py

# manually run scenic_pipeline_downsample_cells.ipynb to downsample the data to four equally sized groups

for rep in {1..4};
do
    cd /scenic_plus_vol/scenic/downsampled_scplus_pipeline_rep${rep}/Snakemake
    export TMPDIR="/scenic_plus_vol/tmp"
    snakemake --cores 90 --rerun-triggers mtime
done