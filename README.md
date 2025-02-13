# immune-cell-census

The goal of this project is to understand natural variation in the immune systems of healthy individuals. 

_The following structures are only proposed at this point_**

### This Repo contains all code for the analyis of the immune cell census project. The structure is as follows:
    . 
    ├── README.md
    ├── Data_Dissemination
    ├── Shared_Resources (for useful install tips and tricks etc)
    ├── Analyses
    |   ├── CoAccessibility
    │   ├── single_cell processing
    │   │   ├── QC_summary
    │   │   ├── archR
    │   │   │   ├── atac_only
    │   │   │   └── joint
    │   │   ├── scanpy
    │   │   │   ├── rna_only
    │   │   │   └── joint
    │   ├── demuxing
    │   ├── QTL Analyses
    │   │   ├── ATAC
    │   │   └── RNA
    │   ├── Mendelian Randomization
    │   ├── testing (misc scripts)
    .   .
    .   .
    .   .
    │   └── Heritability
    ├── Figures
    │   ├── Fig1
    │   ├── Fig2
    │   ├── Fig3
    │   └── Fig4
    └──
    
    
### Data Located on S3 `s3://ucsf-genomics-prod-yelab-icc/` with the following structure:
    . 
    ├── References 
    │   ├── UCSC gene locations
    │   ├── GWAS catalog info
    │   ├── Previous eQTL data
    .   .
    .   .
    .   .
    │   └── CHIP-Peaks
    ├── Production Data
    │   ├──xhlt1 (colabs metis format)
    |   │   ├──data
    |   |   │   ├──single_cell_GEX
    |   |   |   │   ├──processed
    |   |   |   |   │   └──wells (cellranger outs)
    |   |   |   │   └──raw
    |   |   |   |       └──fastqs
    |   |   │   ├──single_cell_ATAC
    |   |   |   │   └──raw
    |   |   |   |       └──fastqs
    │   ├── Genotypes (raw and processed gtypes, will move to xhti) 
    |   |   ├──processed
    |   |   ├   |──cleaned
    |   |   ├   └──imputed
    |   |   └──raw
    │   ├── cellranger outputs (paired down outs for easy syncing)
    │   ├── Freemuxlet
    │   │   ├── pools (inds in each pool + assignments)
    │   │   ├── per well results (XHLT1-POOL-XX#-SC#)
    │   │   ├── per pool results (XHLT1-POOL-XX#)
    |   ├── Processed (junk that should probably be deleted)
    |   ├── Analyses outputs
    │   │   ├── current_objects
    │   |   │   ├── rna
    │   |   │   ├── atac
    │   |   │   └── joint
    │   │   ├── Pseudobulks
    │   │   ├── meta
    │   │   ├── temp_plots (place to stage plots for download while optimizing analyses (make a subdir with your name)
    │   │   │   ├──gracie
    │   │   ├── full genotypes (VCFs + dosage)
    │   │   ├── Processed
    │   │   │   ├── arrowfiles
    │   │   │   ├── Per_wellQC
    │   │   ├── QTL_results
    │   │   │   ├── ATAC
    │   │   │   └── RNA
    │   │   ├── Mendelian Randomization
    .   .
    .   .
    .   .
    │   │   └── H2 files
    ├── Other
    │   ├── Pilot
    │   │   ├── Fastq_Files
    │   │   ├── Cellranger_Outputs
    │   │   └── Freemuxlet
    └── └── 
    
### Protocol for contributing to this body of work
All analyses contributing to this paper will be organized through this repository and performed with AWS with the following procedure. 
 
   1. Each person should work on their own AWS instance and corresponding EBS volume. Follow the guidelines for starting an AWS instance and creating and mounting EBS volumes [here](https://github.com/yelabucsf/AWS/blob/master/AWS_readme.md). **NOTE: To ensure stability of the analyses please either make the script generalizable, or make sure to record to location the sources data on S3.** 

   2. Download this repository to the AWS instance `git clone https://github.com/yelabucsf/immune_cell_census.git`
    
   3. Ensure everything is up to date `git pull`

   4. Create a branch to work on. This should be descriptive of the analysis you are performing and unique to you `git checkout -b my_branch`

   5. Create a new folder in the `Analyses` directory for the analysis you will be performing

   6. Pull data onto your volume. For example, if you are pulling the h5ad files for your analysis please copy the full path as follows. `aws s3 cp --recursive s3://ucsf-genomics-prod-yelab-icc/analysis_outputs/h5ad_files/ /icc/analysis_outputs/h5ad_files/`. This should then be referenced in the code.

   7. Work on your analysis. As you progress save your code and analysis results on your branch using git as follows:
        1. Pull from git to ensure your branch is up to date `git pull`
        2. Stage files for upload `git add myfile`
        3. Commit the changes `git commit -m "Description of changes made"`
        4. Push changes to your branch `git push origin my_branch`
        5. Copy results back to s3 in the folder you created under the `Analyses` directory
        6. **Note: If there are certain files you do not want to ever be added to, such as log files, you can create a `.gitignore` file. This is essentailly just a list of files git will ignore when adding files in bulk. Details can be found [here](https://www.atlassian.com/git/tutorials/saving-changes/gitignore).**

   8. Once you're done with your analysis merge your branch and push to master
        1. After ensuring your branch is up to date following step 7, switch to the master branch `git checkout master`
        2. Make sure your version of master is up to date, `git pull`
        3. Switch back to your branch `git checkout my_branch`
        4. Merge master into your branch `git merge master`
        5. Fix any conflicts git points out. `git status` can be helpful to understand where things may have gone wrong. A more detailed guide can be found [here](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line).
        6. Switch back to master and merge branch to master `git checkout master` then `git merge my_branch`
        7. Push changes to the git server if there are no errors `git push origin master`
        8. If you no longer wish to use your branch, delete it `git branch -d my_branch`
        9. Copy all analysis outputs to s3, `aws s3 sync <analysis dir> s3://ucsf-genomics-prod-yelab-icc/production_analyses/<analysis dir>` Please update the directory structure to include the location of your analysis outputs. **NOTE: all EBS volumes are deleted every two weeks, so it is VERY important to frequently commit your code and sync your analysis outputs to s3/git**

### Think you made a mistake on git? _DON'T PANIC_!
Git tracks changes as you push to it. You can always reset, revert, or return to a previous state in git. A helpful guide can be found [here](https://opensource.com/article/18/6/git-reset-revert-rebase-commands).
