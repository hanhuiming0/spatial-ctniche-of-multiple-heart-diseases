suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
})

MIN_LOGFC <- 0.10
MAX_FDR <- 0.05

files <- c(
  `AMI-IZ` = "5A_DE_AMI_IZ.tsv",
  `ICM-FZ` = "5A_DE_ICM_FZ.tsv",
  HCM = "5A_DE_HCM.tsv"
)
stopifnot(all(file.exists(files)))

de <- lapply(files, read.delim, check.names = FALSE)
up <- lapply(de, function(x) unique(x$gene[x$logFC > MIN_LOGFC & x$FDR < MAX_FDR]))

p <- ggVennDiagram(
  up,
  label_alpha = 0,
  category.names = names(up),
  set_color = c("#D95F5F", "#4C78A8", "#59A14F"),
  set_size = 5,
  label = "count"
) +
  scale_fill_gradient(low = "#F7F7F7", high = "#E68A8A") +
  labs(
    title = "Upregulated genes in ctniche_8 across disease states",
    subtitle = paste0("ctniche_8 vs other niches; log2FC > ", MIN_LOGFC,
                      ", FDR < ", MAX_FDR),
    fill = "Gene count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave("5B_ctniche8_upregulated_DEG_venn.pdf", p, width = 7.5, height = 6.5)

# 输出精确的7个韦恩区域及基因名称。
all_genes <- sort(unique(unlist(up)))
membership <- data.frame(
  gene = all_genes,
  AMI_IZ = all_genes %in% up[["AMI-IZ"]],
  ICM_FZ = all_genes %in% up[["ICM-FZ"]],
  HCM = all_genes %in% up[["HCM"]]
) %>%
  mutate(
    region = case_when(
      AMI_IZ & ICM_FZ & HCM ~ "AMI-IZ & ICM-FZ & HCM",
      AMI_IZ & ICM_FZ ~ "AMI-IZ & ICM-FZ only",
      AMI_IZ & HCM ~ "AMI-IZ & HCM only",
      ICM_FZ & HCM ~ "ICM-FZ & HCM only",
      AMI_IZ ~ "AMI-IZ only",
      ICM_FZ ~ "ICM-FZ only",
      HCM ~ "HCM only"
    )
  )

write.table(membership, "5B_ctniche8_upregulated_DEG_venn_membership.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

counts <- membership %>% count(region, name = "n_genes") %>%
  bind_rows(data.frame(
    region = c("AMI-IZ total", "ICM-FZ total", "HCM total"),
    n_genes = c(length(up[["AMI-IZ"]]), length(up[["ICM-FZ"]]), length(up[["HCM"]]))
  ))
write.table(counts, "5B_ctniche8_upregulated_DEG_venn_counts.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

message("AMI-IZ: ", length(up[["AMI-IZ"]]),
        "; ICM-FZ: ", length(up[["ICM-FZ"]]),
        "; HCM: ", length(up[["HCM"]]),
        "; common: ", length(Reduce(intersect, up)))