

########################
## Script to filter to high quality macs2 peaks from archr project
## peak_dir is the peak_call file in the archr project
## runs code from Ryan Corces
## Author Gracie + Pooja 2025

#conda init
source activate base
conda activate archr


peak_dir="/data/processed/aggregate_phase2_ICC_L1_peaks/PeakCalls/"
out_dir="/data/processed/aggregate_phase2_ICC_L1_peaks/PeakCalls/stringent_peaks"
blacklist_file="/data/processed/meta/reference/hg38-blacklistv2.bed"
prefix="stringent_L1_peaks_spm5e250_rulen12"

mkdir -p ${out_dir}

#generate metadata table
echo "Sample	Group" > ${out_dir}/metadata.txt
ls ${peak_dir}/ReplicateCalls/ | cut -d '-' -f 1 | awk -F'[.]_[.]' '{print $1"._."$2"\t"$1}' >> ${out_dir}/metadata.txt

Rscript createIterativeOverlapPeakSet.R --metadata ${out_dir}/metadata.txt --macs2dir ${peak_dir}/ReplicateCalls/ --outdir ${out_dir}/${prefix} --suffix 'summits.rds' --blacklist ${blacklist_file} --genome hg38 --spm 5 --rule "(n+1)/2" --extend 250
