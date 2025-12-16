# response QTLs in cM

## Environment

This analysis uses the `pyg` environment found in `./env/pyg.txt`.
To install in `./envs/` for example, run the following.
```bash
export ENV_NAME=pyg
mkdir -p ./envs
uv venv ./envs/${ENV_NAME}
source ./envs/${ENV_NAME}/bin/activate
uv pip install -r ./env/pyg.txt
```

The environment can be activated with `./envs/pyg/bin/activate`.

## Data prep

Download the anndata objects for RNA and ATAC.
Filter to include only cM and cM_IL1B cells.
Write the sample IDs present in the anndata.
```bash
./code/cm_il1b/reqtl/subset_adata.py
```

Download the genotypes from s3
```bash
./code/cm_il1b/reqtl/get_geno.sh
```

Filter genotypes for the 333 samples in the anndata file and 
variants with MAF>=5%.
```bash
./code/cm_il1b/reqtl/get_vcfs.sh
```

Pull variant location and ref and alt alleles from the VCF with
```bash
./code/cm_il1b/reqtl/get_snp_loc.sh
```

Pull gene and peak locations into BED file
```bash
python3 ./code/cm_il1b/reqtl/get_feat_loc.py
```

These locations are used for matching peaks/genes and SNPs within cis windows.

## Score cM cells by PC1 as a proxy for IL1B activation

```bash
python3 ./code/cm_il1b/reqtl/score_il1b.py
```

## reQTLs of eQTLs

Get significant eQTLs.
```bash
python3 ./code/cm_il1b/reqtl/filter_eqtls.py
```

Run reQTL analysis.
```bash
python3 ./code/cm_il1b/reqtl/run_reqtl.py
python3 ./code/cm_il1b/reqtl/run_reqtl.il1b.py
```

Plot results
```bash
python3 ./code/cm_il1b/reqtl/plt_reqtl.py
python3 ./code/cm_il1b/reqtl/plt_reqtl.il1b.py
```

## caQTL-specific analysis

This tests where caQTLs are response eQTLs.
The hypothesis is that caQTLs are priming chromatin, and an eQTL can only be detected under certain cell states.

Pull the significant caQTLs in cMs for candidate SNPs.
```bash
python3 ./code/cm_il1b/reqtl/get_cm_caqtls.py
```
It saves the pairs of significant caSNPs and caPeaks.

Now, we test whether the caSNPs are reQTLs for any genes in a cis window.

First, we pull genes with it's TSS +/- 100 kb from the caSNP.
We use 100 kb since calculating reQTLs is expensive.
```bash
python3 ./code/cm_il1b/reqtl/get_cm_casnp_genes.py
```

Run response eQTL analysis.
```bash
python3 ./code/cm_il1b/reqtl/run_reqtl.caqtls.py
```

Run response rcaQTL analysis.
```bash
python3 ./code/cm_il1b/reqtl/run_rcaqtl.caqtls.py
```

## Transfer to S3


```bash
./code/cm_il1b/reqtl/transfer_to_s3.sh
```
