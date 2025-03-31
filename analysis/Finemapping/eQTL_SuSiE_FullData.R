library(nanoparquet)
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(susieR)
library(plinkr)


plink_options <- create_custom_plink_options("/data/Genotype/plink2")


qtl.cutoff <- 5e-4
coloc.windowsize <- 100000
h4.cutoff <- 0.5


this.cell = "B"

this.output.dir = paste0("/data/QTL/RNA/L1_manual/", this.cell, "/")
save.dir = paste0("/data/QTL/RNA/SuSiE/", this.cell, "/")
this.genotype.dir = "/data/Genotype/Cleaned/"

this.cell.pheno = fread(paste0(this.output.dir, "pheno/", this.cell, "_phenos_tensorqtl_sorted.bed"))


for(this.chr in 1:20){
  print(this.chr)
  
  ### This is to make sure that plinkr can process our pgen file where there is no FID
  psam_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.psam")
  df <- data.table::fread(psam_filename, data.table = FALSE)
  names(df) = c("IID", "SEX")
  df$FID = 0
  df = df[, c("FID", "IID", "SEX")]
  # names(df) <- c("FID", names(df)[-1])
  psam_table <- tibble::as_tibble(df)
  psam_table$FID <- as.character(psam_table$FID) # nolint PLINK2 notation
  psam_table$IID <- as.character(psam_table$IID) # nolint PLINK2 notation
  individual_ids <- psam_table$IID # nolint PLINK2 variable name
  
  pvar_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.pvar")
  pvar_table <- plinkr::read_plink2_pvar_file(pvar_filename)
  # pvar_table$fullID = paste0("chr", pvar_table$CHROM, "_", pvar_table$POS, "_", pvar_table$REF, "_", pvar_table$ALT)
  # pvar_table$SNP = ifelse(pvar_table$)
  ### Note: we only map SNPs with rsids 
  snp_names <- pvar_table$ID
  
  pgen_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.pgen")
  pgen_table <- plinkr::read_plink2_pgen_file(
    pgen_filename = pgen_filename,
    names_loci = snp_names,
    names_ind = individual_ids,
    verbose = TRUE
  )
  # attributes(pgen_table)$plinkr_datatype <- "pgen_table"
  
  # this.genotype = read_plink2_pgen_file_from_files(
  #   pgen_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.pgen"),
  #   psam_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.psam"),
  #   pvar_filename = paste0(this.genotype.dir, "chr", this.chr, "_cleaned.pvar"),
  #   verbose = TRUE
  # )
  
  
  
  this.caqtl = read_parquet(paste0(this.output.dir, this.cell, ".cis_qtl_pairs.", this.chr, ".parquet"))
  all.genes = unique(this.caqtl$phenotype_id)
  
  this.cov = fread(paste0(this.output.dir, "cov/", this.cell, "_", this.cell, ".covs"), header = T)
  
  this.individuals = intersect(colnames(this.cov), rownames(pgen_table))
  
  total.nsubject <- length(this.individuals)
  
  to.keep.index = which(colnames(pgen_table) != ".")
  this.geno.matrix = t(pgen_table[this.individuals, to.keep.index])
  rm(pgen_table)
  
  
  
  ### Pheno variance 
  tmp = subset(this.cell.pheno, gene_id %in% all.genes)
  this.pheno.scaled <- as.data.frame(scale(tmp[,5:ncol(tmp)]))
  this.pheno.scaled$gene_id <- tmp$gene_id
  
  
  
  for(this.gene in all.genes){
    print(this.gene)
    # if(file.exists(paste0("ATAC/Pilot/", this.cell, "_", gsub(":", "-", this.peak), "_PIP.csv.gz"))) next
    # print(this.peak)
    this.summary.data = subset(this.caqtl, phenotype_id == this.gene)
    
    
    common.snps = intersect(this.summary.data$variant_id, rownames(this.geno.matrix))
    # common.snps = sort(common.snps)
    
    this.gene.geno = this.geno.matrix[common.snps, ]
    
    this.summary.data = subset(this.summary.data, variant_id %in% common.snps)
    this.summary.data = this.summary.data[match(common.snps, this.summary.data$variant_id), ]
    # this.summary.data$beta_se = this.summary.data$beta / this.summary.data$statistic
    
    # this.peak.geno = subset(this.peak.geno, `CHROM:POS` %in% common.snps)
    
    
    to.remove = which(is.na(this.summary.data$slope_se) | is.na(this.summary.data$slope))
    if(length(to.remove) > 0){
      this.summary.data = this.summary.data[-to.remove,]
      common.snps = common.snps[-to.remove]
    }
    
    ### Require at least 500 SNPs for SuSiE
    if(length(common.snps) < 500) next
    # this.gene.data.pop1 <- subset(pop1.phenos, V1 == this.peak)
    this.gene.data <- subset(this.pheno.scaled, gene_id == this.gene)
    stat.data <- t(as.matrix(this.gene.data[,1:(ncol(this.gene.data)-1)]))
    stat.data <- data.frame(ID = rownames(stat.data),
                            Expression = as.numeric(stat.data))
    
    tmp <- as.data.frame(t(this.cov[,2:ncol(this.cov)]))
    colnames(tmp) <- this.cov$V1
    tmp[,c("Age",paste0("PC", 1:3), paste0("phenPC", 1:10))] <- sapply(tmp[,c("Age",paste0("PC", 1:3), paste0("phenPC", 1:10))], as.numeric)
    tmp$ID <- colnames(this.cov)[2:ncol(this.cov)]
    
    stat.data <- merge(stat.data, tmp, by = c("ID"))
    
    
    # setDT(stat.data)
    this.individuals <- stat.data$ID
    # rownames(stat.data) <- stat.data$ID
    stat.data <- stat.data[,-1]
    stat.data$Sex <- factor(stat.data$Sex)
    stat.data$Site <- factor(stat.data$Site)
    stat.data$GenotypingArray <- factor(stat.data$GenotypingArray)
    
    # stat.data[,c("ID") := NULL]
    stat.data <- stat.data[complete.cases(stat.data),]
    tmp <- lm(Expression ~ ., data = stat.data)
    tmp <- residuals(tmp)
    
    var_y = var(tmp)
    
    
    skip_to_next <- FALSE
    
    tryCatch(       
      {         
        
        # Specifying expression
        # expr = {                      
        #   1 + 1
        #   print("Everything was fine.")
        # },
        
        
        R = cor(t(this.gene.geno[common.snps, ]))
        R[is.na(R)] = 0
        
        fitted_rss1 <- susie_rss(bhat = this.summary.data$slope, 
                                 shat = this.summary.data$slope_se, 
                                 n = total.nsubject, R = R, var_y = var_y, L = 10,
                                 estimate_residual_variance = TRUE)
        
        tmp = data.table(SNP = names(fitted_rss1$pip),
                         PIP = as.numeric(fitted_rss1$pip),
                         CS = 0)
        # tmp$CS[summary(fitted_rss1)$cs$cs] = 1
        if(!is.null(summary(fitted_rss1)$cs)){
          tmp$CS[as.integer(strsplit(summary(fitted_rss1)$cs[["variable"]], ",")[[1]])] = 1  
        }
        
        fwrite(tmp, file = paste0(save.dir, "/", this.cell, "_", this.gene, "_PIP.csv.gz"))
        
      },
      
      # Specifying error message
      error = function(e){          
        print("SuSiE failed")
        skip_to_next <<- TRUE
      }
      
      # warning = function(w){       
      #   print("There was a warning message.")
      # },
      # 
      # finally = {             
      #   print("finally Executed")
      # }
    )
    
    
    if(skip_to_next) { next } 
    
    
  }  
  
  
}





