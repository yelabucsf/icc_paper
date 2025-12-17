import os
import pickle

import numpy as np
from scipy import sparse

import scanpy as sc
import pyranges as pr

from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycistarget.utils import region_names_to_coordinates


import warnings
warnings.filterwarnings("ignore")
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

os.environ['MALLET_MEMORY'] = '256G'


CELL_TYPES = ['cM', 'ncM', 'cDC1', 'cDC2', 'MAIT',  'CD8T_Memory', 'gdT', 
                      'NK_dim', 'NK_bright', 'CD4T_Naive', 'CD8T_Naive', 'Prolif',
                       'CD4T_Memory',  'Treg', 'Plasma', 'B_atypical', 'B_Naive',
                        'B_Memory',  'pDC',  'Progen']
RNA_H5 = "/data_vol/FINAL/SingleCellObjects/rna/MultiVI_posttraining_RNA_Celltypist_VaryingClusteringResolution_manual_annots_filtered.h5ad"
ATAC_H5 = "/data_vol/FINAL/SingleCellObjects/atac/aggregate_phase2_ICC_L2_peaks/aggregate_phase2_ICC_stringent_L2_peaks_counts_mtx.h5"


def main():
    data_dir = "/scenic_plus_vol/scenic/full_dataset_stringent_peaks/100_000_cells"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        os.makedirs(f"{data_dir}/tmp")

    if not os.path.exists(f"{data_dir}/cistopic_obj.pkl"):
        print("Reading in atac adata")
        atac_ad = sc.read_h5ad(ATAC_H5)
        atac_ad.obs.index = atac_ad.obs.index.str.replace("#", ";")  # barcodes have slightly different formatting

        print("Reading in rna adata")
        if not os.path.exists(f"{data_dir}/GEX.h5ad"):
            rna_ad = sc.read_h5ad(RNA_H5)
            rna_ad = rna_ad[rna_ad.obs["L2_manual"].isin(CELL_TYPES)]

            # modify gene names to be ensembl ids
            rna_ad.var["Ensembl ID"] = rna_ad.var.index.values
            rna_ad.var_names = rna_ad.var["name"].values.astype("str")
            rna_ad.var_names_make_unique()
            rna_ad.raw.var["Ensembl ID"] = rna_ad.var["Ensembl ID"].values
            rna_ad.raw.var.index = rna_ad.var_names.values
        else:
            rna_ad = sc.read_h5ad(f"{data_dir}/GEX.h5ad")

        # subset to cells that are present in both rna and atac data
        final_cell_barcodes = np.intersect1d(rna_ad.obs.index.values, atac_ad.obs.index.values)

        print("Subsetting rna adata to cells that are present in both rna and atac data and saving new file")
        rna_ad = rna_ad[rna_ad.obs.index.isin(final_cell_barcodes)]
        if not os.path.exists(f"{data_dir}/GEX.h5ad"):
            sc.write(f"{data_dir}/GEX.h5ad", rna_ad)
        rna_obs = rna_ad.obs
        del rna_ad

        print("Subsetting atac adata to cells that are present in both rna and atac data and adding metadata annotations to atac adata")
        atac_ad = atac_ad[atac_ad.obs.index.isin(final_cell_barcodes)]
        atac_ad.obs = atac_ad.obs.join(rna_obs, how="left") # add metadata annotations to atac adata

        print("Running pycistopic")
        run_pycistopic(data_dir, atac_ad)
    else:
        print("Running pycistopic")
        run_pycistopic(data_dir)

def run_pycistopic(data_dir, atac_ad=None):
    print("Creating cistopic object")
    if not os.path.exists(f"{data_dir}/cistopic_obj.pkl"):
        cistopic_obj = create_cistopic_object(fragment_matrix=sparse.csr_matrix(atac_ad.X.T).astype(np.int32), 
                                        cell_names=atac_ad.obs.index.values,
                                        region_names=atac_ad.var.index.values)
        # Adding cell information
        cell_data = atac_ad.obs
        cell_data.index = atac_ad.obs.index
        cistopic_obj.add_cell_data(cell_data)

        pickle.dump(cistopic_obj,
                    open(f"{data_dir}/cistopic_obj.pkl", "wb"))
    else:
        cistopic_obj = pickle.load(open(f"{data_dir}/cistopic_obj.pkl", "rb"))
    
    if not os.path.exists(data_dir + '/models'):
        os.makedirs(data_dir + '/models')

    print("Running cgs models with mallet")
    if not os.path.exists(f"{data_dir}/models/models_150_iter_LDA.pkl"):
        models=run_cgs_models_mallet(mallet_path="/data/tools/Mallet/bin/mallet",
                            cistopic_obj=cistopic_obj,
                            n_topics=[16],
                            n_cpu=84,
                            n_iter=150,
                            random_state=555,
                            alpha=50,
                            alpha_by_topic=True,
                            eta=0.1,
                            eta_by_topic=False,
                            save_path=data_dir + '/models',
                            tmp_path = data_dir + '/tmp/',
                            reuse_corpus=True)
        
        pickle.dump(models,
                    open(data_dir + '/models/models_150_iter_LDA.pkl', 'wb'))
    else:
        models = pickle.load(open(data_dir + '/models/models_150_iter_LDA.pkl', 'rb'))

    if len(models) > 1:
        model = evaluate_models(models,
                        select_model=16,
                        return_model=True,
                        metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                        plot_metrics=False)
    else:
        model = models[0]
    
    cistopic_obj.add_LDA_model(model)
    # run_umap(cistopic_obj, target  = 'cell', scale=True)
    if "cell" not in cistopic_obj.projections:
        cistopic_obj.projections["cell"] = {}
    if "region" not in cistopic_obj.projections:
        cistopic_obj.projections["region"] = {}
    pickle.dump(cistopic_obj,
                open(f"{data_dir}/cistopic_obj.pkl", 'wb'))
    
    print("Identifying candidate enhancers")
    if not os.path.exists(data_dir + '/candidate_enhancers'):
        os.makedirs(data_dir + '/candidate_enhancers')

    if not os.path.exists(data_dir + '/candidate_enhancers/region_bin_topics_otsu.pkl'):
        region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
        pickle.dump(region_bin_topics_otsu, open(data_dir + '/candidate_enhancers/region_bin_topics_otsu.pkl', 'wb'))
    else:
        region_bin_topics_otsu = pickle.load(open(data_dir + '/candidate_enhancers/region_bin_topics_otsu.pkl', 'rb'))

    if not os.path.exists(data_dir + '/candidate_enhancers/region_bin_topics_top3k.pkl'):
        region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)
        pickle.dump(region_bin_topics_top3k, open(data_dir + '/candidate_enhancers/region_bin_topics_top3k.pkl', 'wb'))
    else:
        region_bin_topics_top3k = pickle.load(open(data_dir + '/candidate_enhancers/region_bin_topics_top3k.pkl', 'rb'))

    # if not os.path.exists(data_dir + '/candidate_enhancers/markers_dict.pkl'):
    #     imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
    #     normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    #     variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
    #     markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='L2_manual', var_features=variable_regions)
    #     pickle.dump(markers_dict, open(data_dir + '/candidate_enhancers/markers_dict.pkl', 'wb'))
    # else:
    #     markers_dict = pickle.load(open(data_dir + '/candidate_enhancers/markers_dict.pkl', 'rb'))
    
    print("Creating region sets")
    region_sets = {}
    region_sets['topics_otsu'] = {}
    region_sets['topics_top_3'] = {}
    # region_sets['DARs'] = {}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    # for DAR in markers_dict.keys():
    #     if len(markers_dict[DAR]) > 0:
    #         regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    #         region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
            
    for k in region_sets.keys():
        if not os.path.exists(f"{data_dir}/candidate_enhancers/{k}"):
            os.makedirs(f"{data_dir}/candidate_enhancers/{k}")
        for topic in region_sets[k].keys():
            region_sets[k][topic].to_bed(f"{data_dir}/candidate_enhancers/{k}/{topic}.bed")


if __name__ == "__main__":
    main()