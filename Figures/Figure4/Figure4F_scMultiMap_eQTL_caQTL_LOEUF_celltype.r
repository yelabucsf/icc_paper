library(data.table)
library(ggplot2)

plot.data <- fread("/data/PaperFigures/Figure4F_scMultiMap_eQTL_caQTL_LOEUF_data.csv.gz")



p = ggplot(plot.data, aes(x = LOEUF, y = propHighPeaks, color = Type)) +
  xlab("LOEUF") +
  ylab("Proportion of Cells with QTL") +
  geom_smooth(
    aes(color = Type),
    method = "glm",
    method.args = list(family = quasibinomial()),
    se = TRUE
  ) +
  stat_cor(
    aes(color = Type),
    method = "pearson",
    # label.sep = "\n",
    label.x.npc = 0.3,
#    label.y.npc = 0.9,
    label.y = seq(0.23, 0.18, length.out = length(unique(plot.data$Type))),
    size = 5
  ) +
  coord_cartesian(ylim = c(0, 0.25)) +
    scale_color_manual(values = c(caqtl.color, eqtl.color)) +
    theme_pubr(base_size = 12, legend = "top") +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) 

ggsave(p, filename = paste0("/data/PaperFigures/Figure4F_celltype_eQTL_caQTL_LOEUF.png"), width = 4.5, height = 4)
ggsave(p, filename = paste0("/data/PaperFigures/Figure4F_celltype_eQTL_caQTL_LOEUF.pdf"), width = 4.5, height = 4)


