library(tidyverse)
library(ChIPseeker)
library(org.Bt.eg.db)
library(clusterProfiler)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)

dmr <- data.table::fread("./dmr.csv")
load("./methylation/all_RNA.Rdata")

module <- data.table::fread("./methylation/wgcna/greenyellow_sub.tpm.tsv")
enhancer <- rbind(
  data.table::fread("./methylation/enhancer/rumen.gene.txt")
)

res <- enhancer[enhancer$GENE_NAME %in% module$V1,]
# res <- enhancer

results <- data.frame()
for (i in 1:nrow(res)) {
  chrom <- str_split(str_split(res[i,3], "[,]", simplify = T), "[:]", simplify = T)[,1]
  locs <- str_split(str_split(res[i,3], "[,]", simplify = T), "[:]", simplify = T)[,2]
  start <- str_split(locs, "[-]", simplify = T)[,1]
  end <- str_split(locs, "[-]", simplify = T)[,2]
  tmp <- data.frame(
    seqnames = chrom,
    start = start,
    end = end,
    gene = res[i, 1]
  )
  results <- rbind(results, tmp)
}

results$seqnames <- sub("chr", "", results$seqnames)
results <- results[order(results$start, results$seqnames),]

res.gr <- makeGRangesFromDataFrame(results, keep.extra.columns = T)
names(dmr)[1] <- "chr"
dmr.gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = T)

final <- GRanges()
for (i in 1:length(res.gr)) {
  tmp <- subsetByOverlaps(dmr.gr, res.gr[i])
  if (length(tmp) > 0) {
    tmp$GENE_NAME <- res.gr[i]$GENE_NAME
    final <- c(final, tmp)
  }
}
final <- as.data.frame(final)
final <- final[!duplicated(final[,1:34]),]
final <- final[grep("adult_rumen", final$hypomethylated_samples),]
data.table::fwrite(as.data.frame(final)[,c(1,2,3,35)], 
                   "feDMR_rumen.bed",
                   quote = F, sep = "\t", row.names = F, col.names = F)
inputs <- final[grep("adult_rumen", final$hypomethylated_samples),
              c("methylation_level_adult_forebrain",
                "methylation_level_adult_muscle",
                "methylation_level_adult_Adipose",
                "methylation_level_adult_rumen",
                "methylation_level_adult_Lung",
                "methylation_level_adult_Spleen")] %>% 
  dplyr::rename("Forebrain" = "methylation_level_adult_forebrain",
         "Muscle" = "methylation_level_adult_muscle",
         "Adipose" = "methylation_level_adult_Adipose",
         "Rumen" = "methylation_level_adult_rumen",
         "Lung" = "methylation_level_adult_Lung",
         "Spleen" = "methylation_level_adult_Spleen") %>% 
  na.omit() %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>% 
  reshape2::melt() %>% 
  mutate(group = "DNA meth.")


inputs2 <- data[unique(res$GENE_NAME),]
inputs2$Forebrain <- apply(inputs2[,which(allTraits$group == "forebrain")],1,
                           mean)
inputs2$Muscle <- apply(inputs2[,which(allTraits$group == "Muscle")],1,
                           mean)
inputs2$Lung <- apply(inputs2[,which(allTraits$group == "Lung")],1,
                           mean)
inputs2$Spleen <- apply(inputs2[,which(allTraits$group == "Spleen")],1,
                           mean)
inputs2$Rumen <- apply(inputs2[,which(allTraits$group == "Rumen")],1,
                           mean)
inputs2$Adipose <- apply(inputs2[,which(allTraits$group == "Adipose")],1,
                           mean)

inputs2 <- inputs2[,c(67, 68, 69, 70, 71, 72)] %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>% 
  reshape2::melt() %>% 
  mutate(group = "Gene expr.")


pdf("greenyellow_EXP.boxplot.pdf", 5, 3)
rbind(inputs, inputs2) %>% 
  ggplot(aes(x = variable, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c(
    "DNA meth." = "#ed2224",
    "Gene expr." = "#344ea2"
  )) +
  labs(x = "", y = "z-score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

