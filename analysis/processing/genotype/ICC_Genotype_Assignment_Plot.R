library(ggpubr)
library(data.table)
library(dplyr)
library(NatParksPalettes)
library(Ternary)


pca.data <- fread("1KGAnchor/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd_finalLiftedHg19..eigenvec")
eigenval <- scan("1KGAnchor/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd_finalLiftedHg19..eigenval")
fam.data <- fread("1KGAnchor/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd_finalLiftedHg19.merge.all_phase3.fam")


pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
percentVar <- round(pve$pve)
onekg.metadata <- fread("../igsr_samples.tsv")




##########################
admixture.fam.data = fread("plinkQC/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd.fam")


##### Let's plot PC distribution for 1KG
pca.data$cohort <- ifelse(pca.data$V2 %in% onekg.metadata$`Sample name`, "1KG", "ICC")


###### just in one plot 
plot.data <- pca.data[,c(1:5, 23)]
names(plot.data) <- c("FID", "IID", "PC1", "PC2", "PC3", "Cohort")
plot.data$refpop <- onekg.metadata$`Superpopulation code`[match(plot.data$IID,
                                                                onekg.metadata$`Sample name`)]
plot.data$refpop <- paste0("1KG ", plot.data$refpop)
plot.data$Population <- ifelse(!grepl("NA", plot.data$refpop), plot.data$refpop, plot.data$Cohort)
plot.data <- plot.data[!is.na(plot.data$Population),]

plot.data = subset(plot.data, Population != "1KG EUR,AFR")


p1 = ggdensity(plot.data, x = "PC1", facet.by = "Population", color = "Population", 
               fill = "Population",
               palette = NatParksPalettes::natparks.pals("Yellowstone", 6)) +
  theme_pubclean()

p2 = ggdensity(plot.data, x = "PC2", facet.by = "Population", color = "Population", 
               fill = "Population",
               palette = NatParksPalettes::natparks.pals("Yellowstone", 6)) +
  theme_pubclean()

p3 = ggdensity(plot.data, x = "PC3", facet.by = "Population", color = "Population", 
               fill = "Population",
               palette = NatParksPalettes::natparks.pals("Yellowstone", 6)) +
  theme_pubclean()

p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T)
ggsave(p, file = "CSGP_GDA_MEGA_AncestryAssignment_1KG_PC_Distribution.png", width = 16.5, height = 6)



### Let's plot 1KG only

plot.data.1kg = subset(plot.data, Population != "ICC")

p1 <- ggscatter(plot.data.1kg, x = "PC1", y = "PC2",
                color = "Population", #palette = park_palette("GeneralGrant", 8),
                palette = "npg",
                # shape = "Cohort", size = "Cohort",
                xlab = paste0("PC1: ", percentVar[1], "% variance"),
                ylab = paste0("PC2: ", percentVar[2], "% variance")) +
  scale_shape_manual(values = c(3, 15, 16)) + 
  scale_size_manual(values = c(1, 1.5, 1.5)) +
  # guides(color=guide_legend(title=paste0("Self-Reported ", this.ancestry))) +
  theme_pubclean() +
  theme(axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13))


p2 <- ggscatter(plot.data.1kg, x = "PC1", y = "PC3",
                color = "Population", palette = "npg",
                # shape = "Cohort", size = "Cohort",
                xlab = paste0("PC1: ", percentVar[1], "% variance"),
                ylab = paste0("PC3: ", percentVar[3], "% variance")) +
  scale_shape_manual(values = c(3, 15, 16)) + 
  scale_size_manual(values = c(1, 1.5, 1.5)) +
  # guides(color=guide_legend(title=paste0("Self-Reported ", this.ancestry))) +
  theme_pubclean() +
  theme(axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13))


p <- ggarrange(p1, p2, ncol = 2, common.legend = T)
ggsave(p, filename = paste0("1KGAnchor/CSGP_GDA01_MEGA2_5_updatedGenotypePCA_1KGOnly", "_anchoredWith1KG.png"),
       width = 16, height = 8)




# ########################################
# ## Remove related individuals
# ########################################
# 

repeated.inds = fread("CGSP_GDA_MEGA_IDmapping_demographics.csv")

to.remove = repeated.inds$IID[which(repeated.inds$RemoveRelated == 1 | repeated.inds$RemoveRepeats == 1)]

plot.data.cleaned = subset(plot.data, !IID %in% to.remove)
plot.data.cleaned = subset(plot.data.cleaned, Population == "ICC")


# ########################################
# Now let's try to resolve admixture
# Four ancestries
# ########################################

### Reassign using 4 ancestryes

tbl=fread("admixture/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd-updated_pruned.4.Q")
tbl$ID = admixture.fam.data$V2

plot.data.cleaned$V1 = tbl$V1[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$V2 = tbl$V2[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$V3 = tbl$V3[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$V4 = tbl$V4[match(plot.data.cleaned$IID, tbl$ID)]


tmp = apply(plot.data.cleaned[, 9:12], 1, which.max)
plot.data.cleaned$Assignment = paste0("V", tmp)

p1 = ggdensity(plot.data.cleaned, x = "PC1", facet.by = "Assignment", color = "Assignment", 
               fill = "Assignment",
               palette = NatParksPalettes::natparks.pals("Volcanoes", 4)) +
  theme_pubclean()


p2 = ggdensity(plot.data.cleaned, x = "PC2", facet.by = "Assignment", color = "Assignment", 
               fill = "Assignment",
               palette = NatParksPalettes::natparks.pals("Volcanoes", 4)) +
  theme_pubclean()


p3 = ggdensity(plot.data.cleaned, x = "PC3", facet.by = "Assignment", color = "Assignment", 
               fill = "Assignment",
               palette = NatParksPalettes::natparks.pals("Volcanoes", 4)) +
  theme_pubclean()

p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T)
ggsave(p, file = "CSGP_GDA_MEGA_AncestryAssignment_Admixture_4Pops_PC_Distribution.png", width = 16.5, height = 4.5)



### Resolve 
tmp = apply(plot.data.cleaned[, 9:12], 1, max)
tmp.min = apply(plot.data.cleaned[, 9:12], 1, min)


### Let's look at distribution of the max score 
p1 = ggdensity(subset(plot.data.cleaned, Assignment == "V1"), x = "V1", fill = NatParksPalettes::natparks.pals("Volcanoes", 4)[1],
               xlab = "EUR Proportion") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey")

# p2 = ggdensity(subset(plot.data.cleaned, Assignment == "V2"), x = "V2", fill = NatParksPalettes::natparks.pals("Volcanoes", 4)[2],
#                xlab = "AFR Proportion") +
#   geom_vline(xintercept = 0.85, linetype = "dashed", color = "grey")

p3 = ggdensity(subset(plot.data.cleaned, Assignment == "V3"), x = "V3", fill = NatParksPalettes::natparks.pals("Volcanoes", 4)[3],
               xlab = "EAS Proportion") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey")

p = ggarrange(p1, p3, ncol = 1, nrow = 2, common.legend = T)
ggsave(p, file = "CGSP_GDA_MEGA_AncestryAssignment_Admixture_Distribution_EUR_EAS_0.9_Proportion_Density.png", width = 4, height = 4)




############################################
## Choose cutoff for ancestry assignment 
## Baesd on distribution of admixture
############################################
this.cutoff = 0.9

tmp = apply(plot.data.cleaned[, 9:12], 1, max)
admixture.resolve = subset(plot.data.cleaned, tmp < this.cutoff)
admixture.confirmed = subset(plot.data.cleaned, tmp >= this.cutoff)


setDT(admixture.resolve)
admixture.resolve[, c("V5") := NULL]


###############################################
### Map to original three-ancestry proportions 
###############################################

tbl=fread("admixture/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd-updated_pruned.3.Q")
tbl$ID = admixture.fam.data$V2

admixture.resolve$EUR = tbl$V1[match(admixture.resolve$IID, tbl$ID)]
admixture.resolve$AFR = tbl$V2[match(admixture.resolve$IID, tbl$ID)]
admixture.resolve$EAS = tbl$V3[match(admixture.resolve$IID, tbl$ID)]


# Configure plotting area
par(mar = rep(0.3, 4))

# Define a colour spectrum
spectrumBins <- 255 # Number of bins to use
mySpectrum <- rowCol <- hcl.colors(spectrumBins, palette = "viridis")

# Cut our reflectance data into categories
binnedReflectance <- cut(admixture.resolve$PC1, spectrumBins)

# Assign each data point a colour from the spectrum
pointCol <- mySpectrum[binnedReflectance]

# # Define a size range
# maxSize <- 2.4 # Size of largest point, in plotting units
# sizes <- dat$grain_size
# pointSize <- sizes * maxSize / max(sizes)

# Initialize the plot
TernaryPlot(atip = "EUR",
            btip = "AFR",
            ctip = "EAS"
)

# Plot the points
# admixture.resolve$AFR = admixture.resolve$V2 + admixture.resolve

TernaryPoints(admixture.resolve[, c("EUR", "AFR", "EAS")],
              # cex = pointSize, # Point size
              col = pointCol,  # Point colour
              pch = 16         # Plotting symbol (16 = filled circle)
)

# CGSP_GDA_MEGA_AncestryAssignment_Admixture_Ternary_0.85.png




######################################
##### Prepare for export 
######################################

### Choose cutoff AMR vs AFR baesd on admixture distribution

afr.cutoff = 0.3
amr.new = admixture.resolve$IID[which(admixture.resolve$AFR <= afr.cutoff)]
afr.new = admixture.resolve$IID[which(admixture.resolve$AFR > afr.cutoff)]
plot.data.cleaned$Assignment[which(plot.data.cleaned$IID %in% amr.new)] = "V4"
plot.data.cleaned$Assignment[which(plot.data.cleaned$IID %in% afr.new)] = "V2"

plot.data.cleaned$Assignment = ifelse(plot.data.cleaned$Assignment == "V1", "EUR",
                                      ifelse(plot.data.cleaned$Assignment == "V2", "AFR",
                                             ifelse(plot.data.cleaned$Assignment == "V3", "EAS", "AMR")))

plot.data.cleaned[, c("V1", "V2", "V3", "V4", "V5", "refpop", "Population") := NULL]


tbl=fread("admixture/CGSP_GDA01_MEGA2_5-postHWE-removeRepeatInd-updated_pruned.4.Q")
tbl$ID = admixture.fam.data$V2

plot.data.cleaned$EUR = tbl$V1[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$AFR = tbl$V2[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$EAS = tbl$V3[match(plot.data.cleaned$IID, tbl$ID)]
plot.data.cleaned$AMR = tbl$V4[match(plot.data.cleaned$IID, tbl$ID)]

plot.data.cleaned$SAS_Admixture = ifelse(plot.data.cleaned$IID %in% sas.pop.inds, 1, 0)
fwrite(plot.data.cleaned, file = "ICC_Genotype_Ancestry_Assignment.csv")


### Export list of individuals for ancestry-specific genotype and QTL analysis
for(this.anc in c("EUR", "AFR", "EAS", "AMR")){
  this.individuals = plot.data.cleaned$IID[which(plot.data.cleaned$Assignment == this.anc)]
  write.table(this.individuals, file = paste0("ICC_", this.anc, "_individuals.txt"),
              row.names = F, col.names = F, sep = "\n", quote = F)
}



######################################
## plot final assignment pie chart
######################################
count.data <- data.frame(table(plot.data.cleaned$Assignment))
names(count.data) <- c( "Population", "Count")
# plot.data$Sex <- ifelse(plot.data$SexLabel == "1", "Male", "Female")

p1 <- ggpie(count.data, x = "Count", label = "Count",
            lab.pos = "in",
            fill = "Population", color = "white",
            palette = NatParksPalettes::natparks.pals("Volcanoes", 4)) +
  theme_pubclean() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(p1, file = "ICC_Ancestry_Assignment_Piechat.png", width = 6, height = 5)













