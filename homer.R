library(tidyverse)
library(dendextend)
library(org.Bt.eg.db)
library(clusterProfiler)



stemCells <- data.table::fread("./methylation/stemCells.tpm.tsv")
tissues <- data.table::fread("./methylation/GSE137943_UMD3.1.1.fpkm.txt")
tissues[,2:ncol(tissues)] <- fpkmToTpm(tissues[,2:ncol(tissues)])
genes <- bitr(stemCells$V1, fromType = "SYMBOL", toType = "ENSEMBL",
              OrgDb = org.Bt.eg.db)
genes <- genes[!duplicated(genes$SYMBOL),]

data <- inner_join(inner_join(genes, stemCells, by=c("SYMBOL"="V1")),
                   tissues, by=c("ENSEMBL"="V1")) %>% 
  remove_rownames() %>% 
  dplyr::select(-c("ENSEMBL")) %>% 
  column_to_rownames(var = "SYMBOL")

allTraits <- data.frame(
  row.names = names(data),
  group = c(
    rep("bEPSCs_AGS", 2),
    rep("bEPSCs_B18", 3),
    rep("bESCs", 15),
    rep("Adipose", 4),
    rep("forebrain", 4),
    rep("Heart", 4),
    rep("Ileum", 4),
    rep("Kidney", 4),
    rep("Muscle", 4),
    rep("Liver", 4),
    rep("Lung", 4),
    "Lymphocyte",
    rep("Mammary", 2),
    rep("Ovary", 2),
    rep("Rumen", 4),
    rep("Spleen", 4),
    rep("Uterus", 2),
    rep("Mammary", 2),
    rep("Lymphocyte", 2)
  )
)

data.log <- log10(data + 1)

inputs <- data.log %>% 
  t() %>% 
  as.data.frame() %>% 
  group_by(allTraits$group) %>% 
  summarise_all(mean)

homer <- data.table::fread("./methylation/homer/all/hypo_kidney/knownResults.txt")
genes <- data.table::fread("./methylation/wgcna/black_sub.tpm.tsv")

intersect(genes$V1, 
          str_to_upper(str_split(homer$`Motif Name`, "[(]", simplify = T)[,1]))

inputs %>% 
  select(c("allTraits$group", "PAX3")) %>% 
  reshape2::melt() %>% 
  rename("sample" = "allTraits$group") %>% 
  ggplot(aes(x = sample, y = value, fill = variable)) +
  geom_col(position="dodge")

x <- inputs %>% 
  remove_rownames() %>% 
  column_to_rownames("allTraits$group") %>% 
  t() %>% 
  as.data.frame() %>% 
  select("bEPSCs_AGS", "bEPSCs_B18", "bESCs", "forebrain",
         "Heart", "Muscle", "Lung", "Adipose", "Uterus", 
         "Rumen", "Kidney", "Ileum",
         "Spleen", "Mammary")
x <- x[c("POU5F1", "NANOG", "TFAP2C", "MYB", "ZIC2", "TEAD2", 
         "SIX3", "SOX10", "NPAS2", "MEF2C",
         "GATA4", "SIX2", "ATF3", "FOXA1", "ERG", "GATA2",  
         "KLF5", "ASCL2", "PRDM1", "PPARA", "HOXD4", "HOXB7", "VDR","HOXB4", 
         "SPDEF", "ELF3", "MEF2B", "IRF1", "IRF2", "IRF8", "BATF", "TCF7"),]
names(x) <- c("hypo_bepscs_ags", "hypo_bepscs_b18", "hypo_bescs_f7",
              "hypo_forebrain", "hypo_heart", "hypo_muscle", "hypo_lung", 
              "hypo_adult_adipose", "hypo_adult_uterus", 
              "hypo_rumen", "hypo_kidney", "hypo_adult_ileum", "hypo_adult_spleen", 
              "hypo_adult_mammary")
for (gene in rownames(x)) {
  for (sample in c("hypo_bepscs_ags", "hypo_bepscs_b18", "hypo_bescs_f7",
                   "hypo_forebrain", "hypo_heart", "hypo_muscle", "hypo_lung", 
                   "hypo_adult_adipose", "hypo_adult_uterus", 
                   "hypo_rumen", "hypo_kidney", "hypo_adult_ileum", "hypo_adult_spleen", 
                   "hypo_adult_mammary")) {
    tmp <- data.table::fread(paste0("./methylation/homer/all/", sample, "/knownResults.txt"))
    if (!gene %in% str_to_upper(str_split(tmp$`Motif Name`, "[(]", simplify = T)[,1])) {
      x[gene, sample] <- NA
    }
  }
}

de.novo <- c("POU5F1", "NANOG", "TFAP2C", "MYB", "ZIC",
             "TEAD2", "SOX10", "GATA4", "SIX2", "FOXA1", 'ERG', 'GATA2',
             "KLF5", "ELF3", "IRF1", "IRF2", "IRF8")
no.de.novo <- c("SIX3", "NPAS2", "MEF2C", "ATF3", "ASCL2", "PRDM1",
                "PPARA", "HOXD4", "HOXB7", "VDR", "HOXB4", "SPDEF",
                "MEF2B", "BATF", "TCF7")
row.anno <- data.frame(
  row.names = rownames(x), 
  class = c("de.novo", "de.novo", "de.novo", "de.novo", "de.novo",
            "de.novo", "non", "de.novo", "non", "non", "de.novo", "de.novo",
            "non", "de.novo", "de.novo", "de.novo", "de.novo","non",  "non",
            "non", "non", "non", "non", "non", "non", "de.novo",
            "non", "de.novo", "de.novo", "de.novo", "non", "non")
)

ann_colors=list(class=c(de.novo='black', non="grey"))
pdf("homer_motif.heatmap.pdf", 8, 4)
pheatmap::pheatmap(t(x), scale = "column", cluster_rows = F, cluster_cols = F,
                   cellwidth = 8, cellheight = 8, border_color = "white",
                   color = colorRampPalette(c("#334fa2","white","#ed2224"))(100),
                   annotation_col = row.anno, annotation_colors = ann_colors)
dev.off()

data <- data.table::fread("./dmr.csv")
data$`#chr` <- paste0("chr", data$`#chr`)
data <- data[data$hypomethylated_samples != "",1:3]
data$name <- paste0("all_", rownames(data))
data$value <- 0
data$strand <- "."

data.table::fwrite(data, "data.bed", quote = F, sep = "\t",
                   row.names = F, col.names = F)
