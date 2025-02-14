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
 
   1. Each person should work on their own AWS instance and corresponding EBS volume. Follow the guidelines for starting an AWS instance and creating and mounting EBS volumes [here](https://github.com/yelabucsf/AWS/blob/master/AWS_readme.md). **NOTE: To ensure stability of the analyses please either make the script generalizable, or make sure to record the source data paths on S3.** 

   2. Download this repository to the AWS instance `git clone https://github.com/yelabucsf/icc_paper.git`
    
   3. Ensure everything is up to date `git pull`

   4. Create a unique branch to work on. This could be your name, or be descriptive of the analysis you are performing `git checkout -b my_branch`

   5. Create a new folder in the `Analyses` directory for the analysis you will be performing

   6. Pull data onto your volume. For example, if you are pulling the h5ad files for your analysis please copy the full path as follows. `aws s3 cp --recursive s3://ucsf-genomics-prod-yelab-icc/analysis_outputs/h5ad_files/ /icc/analysis_outputs/h5ad_files/`. This should then be referenced in the code.

   7. Once you've completed a logical section of your analysis, such as QTL calling, heritability analysis, Figure 2b, etc.
        1. Update the remote github repo with a [pull request](#updating-the-remote-with-pull-requests) or a
           [merge](#updating-the-remote-with-merge).
        2. Copy all analysis outputs to s3, `aws s3 sync <analysis dir> s3://ucsf-genomics-prod-yelab-icc/production_analyses/<analysis dir>`.
           Please update the directory structure to include the location of your analysis outputs.
           **NOTE: all EBS volumes are deleted every two weeks, so it is VERY important to frequently commit your code and sync your analysis outputs to s3/git**
        3. If you no longer wish to use your branch, delete it `git branch -d my_branch`

### Updating the Remote with Pull Requests

TL;DR
1. Create your local working branch (e.g. your name)
2. Add files, make changes, etc. on your working branch
3. Create and merge a PR with main <- your-working-branch
4. Update your local main and working branch

Pull requests are preferred for remote repo updates.

#### Prerequisites

Clone the repo locally
```bash
git clone https://github.com/yelabucsf/icc_paper.git
```

#### 1. Create Your Local Branch

Create your local working branch
```bash
# Create and switch to a new working branch (e.g. your name)
# The -b flag creates a new branch and checks it out in one command
git checkout -b your-working-branch
```

#### 2. Add Files, Make Changes, etc.

Generally keep your local repo up to date and resolve conflicts early

```bash
# Ensure your local main is up to date to avoid future conflicts
git checkout main
git pull origin main
git checkout your-working-branch
git merge main
```

Work on your local working branch
```bash
# Make your changes, then commit
git add .  # Stages all changed files
git commit -m "Descriptive commit message"
# Push and set upstream tracking - the -u flag remembers your branch for future pushes
git push -u origin your-working-branch
```

**Note: If there are certain files you do not want to ever be added to, such as
log files, you can add to the `.gitignore` file. This is essentailly just a list
of files git will ignore when adding files in bulk. Details can be found
[here](https://www.atlassian.com/git/tutorials/saving-changes/gitignore).**

#### 3. Create and Merge a PR

1. **Best Practices**
- One logical change per PR. Makes reviews and reverting easier.
- Split large changes into smaller PRs when possible
  - If you add scripts for Figure 1 and Figure 2, it's easier to create a PR
    for Figure 1 and a separate PR for Figure 2.

2. **Create the PR on GitHub**
- Navigate to the repository on GitHub
- Click "Pull requests" in the top navigation
- Click "New pull request"
- Set base branch (typically `main`) and compare branch (your working branch)
  - The base branch is where you want your changes to end up
  - The compare branch contains your new changes
- Review the changes for accuracy
  - GitHub shows a diff view highlighting additions in green and deletions in red
  - Check that all intended files are included and no unintended files are present
- Click "Create pull request"

3. **PR Components**

Required:
- Title: Clear, concise description of changes. For example, "Add scripts for Sup Fig 1"
- Description: Detailed explanation including what changes were made. You can
  add `@user` to mention and notify someone.

Optional:
- Reviewers: Assign team members if you want them to review your code changes
- Assignees: Usually yourself
- Labels: Categories for the PR (e.g., bug, enhancement), helpful with filtering
- Projects: not very useful
- Milestone: not very useful

4. **Updating Your PR (optional)**

If you forget a file or fix a bug, you can commit these before merging the PR.
The PR automatically updates after new commits to your branch.

```bash
git checkout your-working-branch
# Make additional changes locally
git add .
git commit -m "Fixed a bug in Figure 2a"
git push origin your-working-branch
```

The PR will now show the latest commit and its changes, in addition to
previous commits and changes.

5. **PR Merging Process**

Once you are comfortable with the changes, merge your PR by clicking
"Merge pull request".

**Merge Options Explained**
- Create a merge commit (default)
  - Preserves complete commit history, creating a new commit combining the changes.
- Squash and merge
  - Combines all PR commits into one
  - Useful if you have a ton of commits and just want one for the PR merge.
  - Cleaner history but loses detailed commit context
- Rebase and merge
  - Maintains linear history
  - Moves your commits to the tip of the target branch
  - Useful for keeping a clean project history

Once you confirm the merge and there are no conflicts, that's it! You'll see
the PR was merged and closed. If you go to the main branch, you'll see your
changes. If you click on the commit history, you'll see the merge commit as
well.

6. **Merge Conflicts**

Occurs when the same file and line(s) has been modified in different ways in
your branch and the target branch. Edits to the same file in different lines
don't result in conflicts.
```bash
# Resolve conflicts with main
git checkout your-working-branch
git fetch origin  # Get latest changes
git merge origin/main
# Fix conflicts in your editor, then
git add .
git commit -m "Resolve merge conflicts with main"
git push origin your-working-branch
```
`git status` can be helpful to understand where things may have gone wrong.
A more detailed guide on fixing merge conflicts can be found
[here](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line).

#### 4 Update Your Local Branch

If you don't do this, your working branch won't have changes others have made.
This will likely lead to conflicts down the road.

Update your local main with the merged changes
```bash
# Switch to main and update to get your merged changes
git checkout main
git pull origin main
```

Update your own branch
```bash
# Switch back to working branch
git checkout your-working-branch
git merge main
```

### Updating the Remote with Merge

Although we are trying to use PRs, you could use `pull`, `merge`, and `push`.
```bash
git checkout main
git pull origin main
git checkout your-working-branch
# add your updates and changes...
git commit .
git merge main
# fix any potential conflicts...
git push -u origin your-working-branch
git checkout main
git merge your-working-branch
git push -u origin main
```

### Think you made a mistake on git? _DON'T PANIC_!
Git tracks changes as you push to it. You can always reset, revert, or return to a previous state in git. A helpful guide can be found [here](https://opensource.com/article/18/6/git-reset-revert-rebase-commands).
