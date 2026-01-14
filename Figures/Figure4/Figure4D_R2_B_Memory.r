library(data.table)
library(ggpubr)

plot.data = fread("/data/PaperFigures/Figure4D_R2_data.csv.gz")

p = ggscatter(plot.data, x = "R2_Single", y = "R2",
              #add = "reg.line", conf.int = TRUE,
              xlab = "R2 (Top Peak)", ylab = "R2 (All peaks)",
                color = "caQTL", size = 0.75, #palette = "jco",
                palette = c("eQTL" = "#e41a1c", "caQTL" = "#377eb8", "joint caQTL/eQTL" = "#984ea3"),
            #   title = paste0("L2 ", this.cell, " var explained by ATAC peaks")
            ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    #stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    #scale_color_gradient(low = "blue", high = "red", name = "log10(N Peak + 1)") +
    theme_pubr(base_size = 12, legend = c(0.7, 0.15)) +
    theme(#plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

ggsave(p, file = paste0("/data/PaperFigures/Figure4D_R2_", this.cell, ".png"), width = 4, height = 4)
ggsave(p, file = paste0("/data/PaperFigures/Figure4D_R2_", this.cell, ".pdf"), width = 4, height = 4)



