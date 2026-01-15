### ComplexHeatmap
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(data.table)

plot.data = fread("/data/PaperFigures/Figure3D_coloc_Heatmap_data.csv")

atac.plot.data <- plot.data %>% filter(Assay == "ATAC") %>% pivot_wider(id_cols = Cell, names_from = Trait, values_from = Fraction, values_fill = 0)
atac.mat <- as.matrix(atac.plot.data[,-1])
rownames(atac.mat) <- atac.plot.data$Cell


rna.plot.data <- plot.data %>% filter(Assay == "RNA") %>% pivot_wider(id_cols = Cell, names_from = Trait, values_from = Fraction, values_fill = 0)
rna.mat <- as.matrix(rna.plot.data[,-1])
rownames(rna.mat) <- rna.plot.data$Cell   


ht.atac <- Heatmap(
  atac.mat,
  name = "Coloc Fraction\n(ATAC)",
  col = colorRamp2(c(0, 0.2), c("#eff3ff", "#084594")),
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
ht.atac.drawn <- draw(ht.atac)

row.order <- row_order(ht.atac.drawn)
col.order <- column_order(ht.atac.drawn)

ht.rna <- Heatmap(
  rna.mat,
  name = "Coloc Fraction\n(RNA)",
  col = colorRamp2(c(0, 0.2), c("#fee5d9", "#99000d")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = row.order,
  column_order = col.order
)

ht.combined <- ht.atac %v% ht.rna


pdf(
  "/data/PaperFigures/Figure3D_coloc_Heatmap.pdf",
  width = 5,
  height = 8
)

draw(ht.combined, heatmap_legend_side = "right")

dev.off()

