#!/bin/bash

mamba activate scenicplus
python run_pycistopic.py

# manually run scenic_pipeline_downsample_cells.ipynb to downsample the data to 100,000 cells

cd /scenic_plus_vol/scenic/100_000_cell_scplus_pipeline//Snakemake
export TMPDIR="/scenic_plus_vol/tmp"
snakemake --cores 128