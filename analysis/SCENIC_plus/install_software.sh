# create new env and install scenicplus
mamba create --name scenicplus python=3.11.8 -y
mamba activate scenicplus
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install .

# install mallet
apt-get install -y --no-install-recommends ant openjdk-11-jdk
git clone --depth=1 https://github.com/mimno/Mallet.git
cd Mallet
ant

# download motif database
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather
wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl
# download TF names
wget http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt