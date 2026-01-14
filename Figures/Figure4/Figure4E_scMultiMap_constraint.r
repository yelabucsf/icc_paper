library(data.table)
library(ggpubr)

this.cell = "B_memory"
stat.data <- fread("/data/PaperFigures/Figure4E_scMultiMap_constraint_data.csv.gz")

p = ggplot(stat.data, aes(x = LOEUF_decile, y = nPeaks)) +
    xlab("LOEUF bin") +
    ylab("N Significant Peaks") +
    #scale_x_continuous(breaks = unique.loeuf) +
    geom_boxplot(aes(group = LOEUF_decile), color = this.color, outliers = FALSE) +
    #geom_jitter(color = this.light.color) +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(p.value.label), sep = "~")),
                    #label.y = 3, 
                    coef.digits = 2,
                    parse = TRUE,
                    #rr.digits = 2,
                    p.digits = 3) +
  theme_pubr(base_size =  12) +
#   stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3.5) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
ggsave(p, filename = paste0("/data/PaperFigures/Figure4E_scMultiMap_", this.cell, "_LOEUF_nSigPeaks_boxplot.png"),
       width = 4.5, height = 3)
ggsave(p, filename = paste0("/data/PaperFigures/Figure4E_scMultiMap_", this.cell, "_LOEUF_nSigPeaks_boxplot.pdf"),
       width = 4.5, height = 3)