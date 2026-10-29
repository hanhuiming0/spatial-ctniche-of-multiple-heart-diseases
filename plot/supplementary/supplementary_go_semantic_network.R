suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(cluster)
})

INFILE <- if (file.exists("enrichG.8.txt")) "enrichG.8.txt" else file.path("..", "enrichG.8.txt")
OUT_PREFIX <- "ctniche8_GOBP"
FDR_CUTOFF <- 0.05
NETWORK_SIM_CUTOFF <- 0.70
K_CANDIDATES <- 5:14
N_LABEL_PER_GROUP <- 2

x <- read.delim(INFILE, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(!is.na(p.adjust), p.adjust < FDR_CUTOFF) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    GeneRatio_value = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio)),
    minus_log10_FDR = -log10(pmax(p.adjust, 1e-300)),
    Description_short = ifelse(
      nchar(Description) > 58,
      paste0(substr(Description, 1, 55), "..."),
      Description
    )
  )

if (nrow(x) < 5) stop("<U+663E><U+8457>GO<U+6761><U+76EE><U+4E0D><U+8DB3>5<U+4E2A>")

# 1. <U+57FA><U+4E8E>GO<U+672C><U+4F53><U+7ED3><U+6784><U+8BA1><U+7B97>Wang<U+8BED><U+4E49><U+76F8><U+4F3C><U+6027>
sem_data <- godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE)
sim <- mgoSim(x$ID, x$ID, semData = sem_data, measure = "Wang", combine = NULL)
sim[is.na(sim)] <- 0
diag(sim) <- 1
sim <- (sim + t(sim)) / 2

dist_mat <- 1 - sim
dist_mat[dist_mat < 0] <- 0
d <- as.dist(dist_mat)
hc <- hclust(d, method = "average")

# <U+5728><U+5408><U+7406><U+8303><U+56F4><U+5185><U+7528><U+5E73><U+5747>silhouette<U+9009><U+62E9><U+5927><U+7EC4><U+6570><U+91CF><U+FF0C><U+800C><U+4E0D><U+662F><U+56FA><U+5B9A><U+624B><U+5DE5><U+5206><U+7EC4><U+3002>
k_use <- K_CANDIDATES[K_CANDIDATES < nrow(x)]
sil <- sapply(k_use, function(k) {
  cl <- cutree(hc, k = k)
  mean(cluster::silhouette(cl, d)[, "sil_width"])
})
best_k <- if (length(sil) > 0 && any(is.finite(sil))) k_use[which.max(sil)] else 8L
x$Functional_group_id <- cutree(hc, k = best_k)[x$ID]

# <U+6BCF><U+7EC4><U+4EE3><U+8868><U+6761><U+76EE><U+FF1A><U+517C><U+987E><U+7EC4><U+5185><U+4E2D><U+5FC3><U+6027><U+3001><U+663E><U+8457><U+6027><U+548C><U+5BCC><U+96C6><U+57FA><U+56E0><U+6570><U+3002>
representatives <- bind_rows(lapply(sort(unique(x$Functional_group_id)), function(g) {
  idx <- which(x$Functional_group_id == g)
  centrality <- if (length(idx) == 1) 1 else rowMeans(sim[idx, idx, drop = FALSE])
  zc <- as.numeric(scale(centrality))
  zc[!is.finite(zc)] <- 0
  zs <- as.numeric(scale(x$minus_log10_FDR[idx]))
  zs[!is.finite(zs)] <- 0
  zn <- as.numeric(scale(log1p(x$Count[idx])))
  zn[!is.finite(zn)] <- 0
  score <- zc + 0.35 * zs + 0.20 * zn
  x[idx[which.max(score)], c("Functional_group_id", "ID", "Description")] %>%
    dplyr::rename(Representative_GO = ID, Functional_group = Description)
}))

x <- x %>%
  left_join(representatives, by = "Functional_group_id") %>%
  mutate(
    Functional_group_label = paste0("G", Functional_group_id, ": ", Functional_group),
    Functional_group_label = factor(
      Functional_group_label,
      levels = representatives %>% arrange(Functional_group_id) %>%
        transmute(z = paste0("G", Functional_group_id, ": ", Functional_group)) %>% pull(z)
    ),
    # <U+5C06><U+81EA><U+52A8><U+8BED><U+4E49><U+7C07><U+8FDB><U+4E00><U+6B65><U+5F52><U+5E76><U+4E3A><U+8BBA><U+6587><U+4E3B><U+56FE><U+53EF><U+89E3><U+91CA><U+7684>5<U+4E2A><U+751F><U+7269><U+5B66><U+5927><U+7C7B><U+3002>
    # <U+8BE5><U+6620><U+5C04><U+9488><U+5BF9><U+5F53><U+524D>enrichG.8.txt<U+FF1B><U+82E5><U+8F93><U+5165><U+7ED3><U+679C><U+6539><U+53D8><U+FF0C><U+5E94><U+91CD><U+65B0><U+6838><U+5BF9><U+4EE3><U+8868><U+6761><U+76EE><U+3002>
    Supergroup = case_when(
      Functional_group_id %in% c(1, 2, 3, 4) ~ "ECM, collagen and mesenchymal remodelling",
      Functional_group_id %in% c(5, 7, 8) ~ "TGF-beta and receptor signalling",
      Functional_group_id == 9 ~ "Hypoxia and oxygen response",
      Functional_group_id %in% c(6, 10, 12, 13) ~ "Immune and inflammatory response",
      Functional_group_id %in% c(11, 14) ~ "Protease regulation and extracellular secretion",
      TRUE ~ "Other"
    )
  )

write.table(
  x %>% arrange(Functional_group_id, p.adjust),
  paste0(OUT_PREFIX, "_all_terms_with_groups.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 2. <U+5927><U+529F><U+80FD><U+7EC4><U+6982><U+89C8><U+56FE>
group_summary <- x %>%
  group_by(Functional_group_id, Functional_group_label, Functional_group) %>%
  summarise(
    n_terms = n(),
    n_unique_genes = n_distinct(unlist(strsplit(paste(geneID, collapse = "/"), "/"))),
    best_FDR = min(p.adjust),
    median_fold_enrichment = median(FoldEnrichment),
    .groups = "drop"
  )

write.table(group_summary, paste0(OUT_PREFIX, "_group_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

supergroup_summary <- x %>%
  group_by(Supergroup) %>%
  summarise(
    n_terms = n(),
    n_unique_genes = n_distinct(unlist(strsplit(paste(geneID, collapse = "/"), "/"))),
    best_FDR = min(p.adjust),
    median_fold_enrichment = median(FoldEnrichment),
    .groups = "drop"
  )
write.table(supergroup_summary, paste0(OUT_PREFIX, "_supergroup_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p_super <- ggplot(
  supergroup_summary,
  aes(
    median_fold_enrichment,
    reorder(Supergroup, median_fold_enrichment)
  )
) +
  geom_segment(aes(x = 0, xend = median_fold_enrichment, yend = Supergroup),
    colour = "grey80", linewidth = 0.5
  ) +
  geom_point(aes(size = n_terms, colour = -log10(best_FDR))) +
  scale_colour_gradient(low = "#4C78A8", high = "#B2182B") +
  scale_size_continuous(range = c(5, 13)) +
  labs(
    title = "ctniche_8 major biological programmes",
    x = "Median fold enrichment", y = NULL,
    colour = expression(-log[10](best ~ FDR)), size = "GO terms"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
ggsave(paste0(OUT_PREFIX, "_supergroup_overview.pdf"), p_super, width = 10.5, height = 4.8)

p_group <- ggplot(
  group_summary,
  aes(x = median_fold_enrichment, y = reorder(Functional_group_label, median_fold_enrichment))
) +
  geom_segment(aes(x = 0, xend = median_fold_enrichment, yend = Functional_group_label),
    colour = "grey80", linewidth = 0.5
  ) +
  geom_point(aes(size = n_terms, colour = -log10(best_FDR))) +
  scale_colour_gradient(low = "#4C78A8", high = "#B2182B") +
  scale_size_continuous(range = c(4, 12)) +
  labs(
    title = "ctniche_8 GOBP functional groups",
    x = "Median fold enrichment", y = NULL,
    colour = expression(-log[10](best ~ FDR)), size = "GO terms"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))

ggsave(paste0(OUT_PREFIX, "_group_overview.pdf"), p_group, width = 10.5, height = 6.5)

# 3. <U+6240><U+6709>224<U+6761><U+529F><U+80FD><U+7684><U+5206><U+9875><U+56FE><U+FF1A><U+6BCF><U+4E2A><U+8BED><U+4E49><U+5927><U+7EC4><U+4E00><U+9875><U+FF0C><U+7EDD><U+4E0D><U+4E22><U+5F03><U+6761><U+76EE>
pdf(paste0(OUT_PREFIX, "_all_terms_by_group.pdf"),
  width = 12, height = 9,
  onefile = TRUE
)
for (g in sort(unique(x$Functional_group_id))) {
  z <- x %>%
    filter(Functional_group_id == g) %>%
    arrange(FoldEnrichment) %>%
    mutate(Description_page = factor(Description, levels = unique(Description)))

  p <- ggplot(z, aes(FoldEnrichment, Description_page)) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, yend = Description_page),
      colour = "grey85", linewidth = 0.4
    ) +
    geom_point(aes(size = Count, colour = minus_log10_FDR)) +
    scale_colour_gradient(low = "#4575B4", high = "#A50026") +
    scale_size_continuous(range = c(2.5, 9)) +
    labs(
      title = unique(as.character(z$Functional_group_label)),
      subtitle = paste0(nrow(z), " enriched GOBP terms"),
      x = "Fold enrichment", y = NULL,
      colour = expression(-log[10](FDR)), size = "Gene count"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = ifelse(nrow(z) > 35, 6.5, 8))
    )
  print(p)
}
dev.off()

# 4. GO<U+8BED><U+4E49><U+7F51><U+7EDC><U+FF1A><U+8282><U+70B9>=GO<U+6761><U+76EE><U+FF1B><U+8FB9>=Wang<U+8BED><U+4E49><U+76F8><U+4F3C><U+6027>>=<U+9608><U+503C>
# <U+4E3A><U+907F><U+514D>224<U+4E2A><U+6807><U+7B7E><U+4E92><U+76F8><U+8986><U+76D6><U+FF0C><U+4EC5><U+6807><U+6CE8><U+6BCF><U+7EC4><U+6700><U+4E2D><U+5FC3>/<U+663E><U+8457><U+7684><U+524D>2<U+6761><U+FF0C><U+6240><U+6709><U+8282><U+70B9><U+4ECD><U+4FDD><U+7559><U+3002>
edge_idx <- which(upper.tri(sim) & sim >= NETWORK_SIM_CUTOFF, arr.ind = TRUE)
edges <- data.frame(
  from = rownames(sim)[edge_idx[, 1]],
  to = colnames(sim)[edge_idx[, 2]],
  similarity = sim[edge_idx]
)
nodes <- x %>%
  transmute(
    name = ID,
    Description,
    Functional_group_id = factor(Functional_group_id),
    Functional_group_label,
    Supergroup,
    Count,
    minus_log10_FDR,
    label_score = minus_log10_FDR + log1p(Count)
  )

label_ids <- nodes %>%
  group_by(Functional_group_id) %>%
  slice_max(label_score, n = N_LABEL_PER_GROUP, with_ties = FALSE) %>%
  pull(name)
nodes$plot_label <- ifelse(nodes$name %in% label_ids, nodes$Description, "")

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
set.seed(20260826)
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = similarity), colour = "grey55", show.legend = FALSE) +
  geom_node_point(aes(
    size = Count, colour = Supergroup,
    shape = minus_log10_FDR > median(minus_log10_FDR)
  ), alpha = 0.88) +
  geom_node_text(aes(label = plot_label),
    repel = TRUE, size = 2.6,
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2
  ) +
  scale_edge_alpha(range = c(0.08, 0.50)) +
  scale_size_continuous(range = c(2.2, 8), name = "Gene count") +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    breaks = c(FALSE, TRUE),
    labels = c(
      `FALSE` = "Lower relative significance (-log10 FDR <= median)",
      `TRUE` = "Higher relative significance (-log10 FDR > median)"
    ),
    name = "Node fill"
  ) +
  labs(
    title = "ctniche_8 GOBP semantic similarity network",
    subtitle = paste0("224 significant terms; edge: Wang similarity >= ", NETWORK_SIM_CUTOFF),
    colour = "Major programme"
  ) +
  theme_graph(base_family = "sans") +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(paste0(OUT_PREFIX, "_semantic_network.pdf"), p_network,
  width = 14, height = 11, limitsize = FALSE
)

# <U+4FDD><U+5B58><U+8BED><U+4E49><U+76F8><U+4F3C><U+6027><U+548C><U+7F51><U+7EDC><U+8FB9><U+FF0C><U+4FBF><U+4E8E><U+8C03><U+6574><U+9608><U+503C><U+800C><U+4E0D><U+91CD><U+590D><U+8BA1><U+7B97><U+3002>
saveRDS(sim, paste0(OUT_PREFIX, "_semantic_similarity_matrix.RDS"))
write.table(edges, paste0(OUT_PREFIX, "_semantic_network_edges.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message(
  "GO<U+5C55><U+793A><U+5B8C><U+6210><U+FF1A>", nrow(x), "<U+6761><U+FF1B><U+81EA><U+52A8><U+8BED><U+4E49><U+5927><U+7EC4><U+FF1A>", best_k,
  "<U+FF1B><U+7F51><U+7EDC><U+8FB9><U+FF1A>", nrow(edges)
)
