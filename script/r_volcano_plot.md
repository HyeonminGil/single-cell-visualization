### cf) https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
```R
## Volcano plot =====
my_res <- res %>%
    as.data.frame() %>%
    filter(!is.na(padj))
dim(my_res)

DEG_res <- read.csv(glue("{dirout}DEG_list.csv"), header = TRUE, row.names = 1)
dim(DEG_res)
head(DEG_res)

vplot_data <- res %>% 
  as.data.frame() %>%
  mutate(
    Expression = case_when(res$padj < 0.05 ~ "Significant",
                           TRUE ~ "Not Significant"),
    logadjp = -log10(res$padj)
    )
head(vplot_data)

vplot_genes <- c(
    "Ctsl", "Atp6v0d2", "Fabp4", "Gpnmb", "Epb41l3", "Fabp5", "Pld3", "Htra1", "Trem2",
    "Mgl2", "Tnfsf9", "Plac8", "Il1b", "Ccr2", "Cx3cr1", "Lifr", "Kmo", "Areg", "Nr4a1", "Ehd1", "Nfkbia")

top_genes <- vplot_data %>%
    filter(rownames(vplot_data) %in% vplot_genes)
head(top_genes)

vplot_data$label <- ""
idx_to_label <- match(vplot_genes, rownames(vplot_data))
vplot_data$label[idx_to_label] <- rownames(vplot_data)[idx_to_label]

vplot_data <- mutate(
    vplot_data, gene_color = case_when(log2FoldChange > 0 ~ "Blue", TRUE ~ "Red")
)

library(ggrepel)
p <- ggplot(vplot_data, aes(log2FoldChange, logadjp, label = label, color = gene_color)) +
  theme_bw() +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("Log2 Fold Change")) + 
  ylab(expression("Adjust P value")) +
  xlim(c(-10, 12)) +
  scale_color_manual(values = c("blue", "black", "red", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel() +
  theme(text = element_text(size = 12),
    axis.text = element_text(family = "Arial", color = "black", size = 12),
    legend.text = element_text(family = "Arial", color = "black", size = 10),
    legend.title = element_text(family = "Arial", color = "black", size = 10))
ggsave(glue("{dirout}3-2_volcano_plot.jpg"), p, width = 6, height = 5)

```
