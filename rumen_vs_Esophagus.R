library(tidyverse)
library(edgeR)

data <- data.table::fread("./methylation/exp_data/all_rumen.exp.counts.tsv") %>% 
  remove_rownames() %>% 
  column_to_rownames("Geneid")

fpkm <- as.data.frame(scuttle::calculateFPKM(data[,2:ncol(data)], data$Length))


group <- c(rep("other", 2), rep("adult_Rumen", 5),
           "other")

design_pheno <- model.matrix(~ 0 + group)
y <- DGEList(counts = data[,2:ncol(data)], genes = rownames(data), group = group)
print(y)

keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design_pheno)

res <- exactTest(y, pair = c("adult_Rumen", "other"))
diff <- topTags(res, n = Inf)$table

EnhancedVolcano::EnhancedVolcano(
  toptable = diff,
  lab = diff$genes,
  x = "logFC", y = "FDR", pCutoff = 0.05, FCcutoff = 1,
  subtitle = "", caption = ""
)

pval <- 0.05
logfc <- 2

diffSig <- diff[(diff$FDR < pval & (diff$logFC > logfc | diff$logFC < -logfc)), ]
write.csv(diffSig, file = "diffSig.csv", quote = FALSE)
diffUp <- diff[(diff$FDR < pval & (diff$logFC > logfc)), ] # for go, kegg
write.csv(diffUp, file = "diffUp.csv", quote = FALSE)
diffDown <- diff[(diff$FDR < pval & (diff$logFC < -logfc)), ] # for go, kegg
write.csv(diffDown, file = "diffDown.csv", quote = FALSE)


fpkm <- data.table::fread("./methylation/exp_data/all_rumen.exp.fpkm.tsv")

inputs <- rbind(
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% data.table::fread("./Esophagus/diffUp.csv")$genes), # Esophagus special
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% data.table::fread("./adult_Rumen/diffDown.csv")$genes), # Esophagus and fetal Rumen special
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% data.table::fread("./fetal_Rumen/diffUp.csv")$genes), # fetal Rumen special
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% data.table::fread("./fetal_Rumen/diffDown.csv")$genes), # fetal Rumen special
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% data.table::fread("./adult_Rumen/diffUp.csv")$genes) # adult Rumen special

)

inputs <- rbind(
  fpkm %>%
    as.data.frame() %>%
    dplyr::filter(V1 %in% diffUp$genes) # Esophagus special
  
)

inputs <- inputs[,
                 c("Esophagus", "fetal_Rumen_1", "fetal_Rumen_2",
                   "adult_Rumen_1", "adult_Rumen_2", "adult_Rumen_3",
                   "adult_Rumen_4", "adult_Rumen_5")]

pdf("heatmap.pdf", 7, 5)
pheatmap::pheatmap(inputs, show_rownames = F, cluster_cols = F, scale = "row",
                   cluster_rows = T, clustering_method = "ward.D2",
                   # gaps_row = c(159, 159 + 728, 159 + 728 + 1106),
                   angle_col = 45, cellwidth = 15,
                   cutree_rows = 6)
dev.off()

row_cluster <- cutree(p$tree_row, k = 6)
newOrder <- inputs[p$tree_row$order,]
inputs <- inputs[which(row_cluster %in% c(4)),]
pheatmap::pheatmap(inputs, show_rownames = F, cluster_cols = F, scale = "row",
                   cluster_rows = T, clustering_method = "ward.D2",
                   # gaps_row = c(159, 159 + 728, 159 + 728 + 1106),
                   angle_col = 45, cellwidth = 15)

res.pca <- prcomp(fpkm[, 2:9], scale = T, retx = T, rank = 3)

inputs <- as.data.frame(res.pca$rotation) %>%
  rownames_to_column(var = "sample")

inputs$group <- c(
  rep("fetal_Rumen", 2), rep("adult_Rumen", 5), "Esophagus"
)

pdf("pca.pdf", 5, 5)
ggplot(inputs, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 5) +
  # ggrepel::geom_text_repel(
  #   data = inputs,
  #   aes(label = sample),
  #   size = 3,
  #   max.overlaps = 200,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")
  # ) +
  cowplot::theme_cowplot()
dev.off()

res.pca$x

# chr3: 16751000-18773000
EDC_genes <- c(
  'S100A1',
  'S100A13',
  'S100A14',
  'S100A16',
  'S100A2',
  'S100A3',
  'S100A4',
  'S100A5',
  'S100A7',
  'S100A8',
  'S100A12',
  'S100A9',
  'PGLYRP4',
  'PGLYRP3',
  'LOR',
  'PRR9',
  'LELP1',
  'SPRR4',
  'IVL',
  'SMCP',
  # 'LCE1B',
  'KPRP',
  'LCE3C',
  'CRCT1',
  'CRNN',
  'RPTN',
  'TCHH',
  'TCHHL1',
  'S100A10',
  'S100A11',
  'S100A10',
  "TGM1",
  'TGM2',
  "TGM3",
  "TGM5",
  "KRT8",
  "KRT18",
  "KRT36",
  "KLK10",
  "KLK12",
  "SPINK5",
  "DUOX1",
  "DUOX2",
  "DUOXA1",
  "DUOXA2",
  "PIP"
)


inputs <- fpkm %>% 
  dplyr::filter(V1 %in% EDC_genes) %>% 
  remove_rownames() %>% 
  column_to_rownames("V1") %>% 
  dplyr::select(c("Esophagus"), everything())

inputs <- log10(inputs + 1)

pdf("EDC.heatmap.pdf", 5, 7)
pheatmap::pheatmap(inputs, cluster_rows = T, cluster_cols = F,
                   clustering_method = "ward.D2",
                   angle_col = 45, border_color = NA, cellwidth = 20, cellheight = 15)
dev.off()


data <- data.table::fread("./fetal_Rumen/diffUp.csv")
gost_fetal <- clusterProfiler::enrichGO(data$genes,
                                        OrgDb = org.Bt.eg.db,
                                        keyType = "SYMBOL",
                                        ont = "BP")
View(as.data.frame(gost_fetal))
# gost_fetal <- gprofiler2::gost(data$genes, organism = "btaurus")$result[,1:13]
# gost_adult <- gprofiler2::gost(data$genes, organism = "btaurus")$result[,1:13]

load("./methylation/exp_data/all_RNA.Rdata")
