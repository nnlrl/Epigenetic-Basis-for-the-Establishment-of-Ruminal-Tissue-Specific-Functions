library(tidyverse)
library(pheatmap)
library(ChIPseeker)
library(org.Bt.eg.db)
library(ComplexHeatmap)
library(clusterProfiler)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
setwd("./methylation/enhancer/")

res <- list()
all <- GRanges()

data <- data.table::fread("Rumen_1.txt")
data.gr <- makeGRangesFromDataFrame(data,
                                    seqnames.field = "CHROM",
                                    start.field = "START",
                                    end.field = "STOP")

res[["Rumen"]] <- data.gr
all <- union(all, data.gr)
for (sample in c("Spleen", "Muscle", "Adipose", "Lung",
                 "Cortex")) {
  
  file_1 <- paste0(sample, "_1.txt")
  file_2 <- paste0(sample, "_2.txt")
  
  data_1 <- data.table::fread(file_1)
  data.gr.1 <- makeGRangesFromDataFrame(data_1,
                                        seqnames.field = "CHROM",
                                        start.field = "START",
                                        end.field = "STOP")
  data_2 <- data.table::fread(file_2)
  data.gr.2 <- makeGRangesFromDataFrame(data_2,
                                        seqnames.field = "CHROM",
                                        start.field = "START",
                                        end.field = "STOP")
  data.gr <- union(data.gr.1, data.gr.2)
  
  res[[sample]] <- data.gr
  all <- union(all, data.gr)
}

data.table::fwrite(as.data.frame(all), "all_enhancer.tsv", sep = "\t", quote = F,
                   row.names = F)
data.table::fwrite(as.data.frame(all)[,1:3], "all_enhancer.bed", sep = "\t", quote = F,
                   row.names = F, col.names = F)

genes <- lapply(res, function(i) {
  seq2gene(i, c(-1000, 3000), 5000, TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene)
})

annopeak <- annotatePeak(all, flankDistance = 5000,
                         TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene)
annopeak.df <- as.data.frame(annopeak)

pdf("superEnhancer.pie.pdf", 7, 7)
plotAnnoPie(annopeak)
dev.off()
plotDistToTSS(annopeak)

genes <- bitr(annopeak.df$geneId, fromType = "ENTREZID", toType = "SYMBOL",
              OrgDb = org.Bt.eg.db)

cc <- compareCluster(geneClusters = genes, 
                     fun="enrichGO", OrgDb = org.Bt.eg.db)               
dotplot(cc, showCategory=10)

length(setdiff(genes[["Rumen"]],
               unlist(genes[2:6])))

genes <- seq2gene(all, c(-1000, 3000), 5000, TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene)
go <- enrichGO(genes, OrgDb = org.Bt.eg.db, ont = "ALL")
dotplot(go)
kegg <- enrichKEGG(genes, organism = "bta")
dotplot(kegg)


# overlap with large hypo CG-DMRs
enhancer <- makeGRangesFromDataFrame(data.table::fread("../enhancer/all_enhancer.bed"),
                                     seqnames.field = "V1",
                                     start.field = "V2",
                                     end.field = "V3")

files <- list.files(pattern = "*.hypo.2kb.bed")


res <- data.frame()
all <- GRanges()
overlap <- data.frame()
for (file in files) {
  data <- data.table::fread(file)
  data$V1 <- paste0("chr", data$V1)
  res <- rbind(res, data)
  data.gr <- makeGRangesFromDataFrame(
    data,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
  )
  overlap <- rbind(overlap, c(file, sum(countOverlaps(data.gr, enhancer)),
                              length(data.gr)))
  
  
  all <- union(all, data.gr)
}
names(overlap) <- c("sample", "V1", "V2")

overlap$percent <- as.numeric(overlap$V1) / as.numeric(overlap$V2)

inputs <- overlap[grep("([Aa]dipose|[Ll]ung|[Mm]uscle|[Rr]umen|Spleen|[Bb]rain)",
                       overlap$sample),]

inputs$group <- c("Adipose", "Forebrain", "Lung", "Muscle",
                  "Rumen", "Spleen", "Forebrain", "Hindbrain",
                  "Lung", "Muscle", "Rumen")
inputs$V1 <- as.numeric(inputs$V1)
inputs$V2 <- as.numeric(inputs$V2)

pdf("enhancer_intersect_largeDMR.barplot.pdf", 6, 4)
inputs %>% 
  select(c("V1", "V2", "group")) %>% 
  group_by(group) %>% 
  summarise_all(sum) %>%
  filter(group %in% c("Adipose", "Spleen", "Lung", 
                      "Muscle", "Rumen", "Forebrain")) %>% 
  rename("Intersection of super-enhancer" = "V2",
         "Other" = "V1") %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = factor(group,
                        levels = c("Adipose",
                                   "Spleen",
                                   "Forebrain",
                                   "Lung",
                                   "Muscle",
                                   "Rumen")), y = value, fill = variable)) +
  geom_col(width = 0.5) +
  # geom_text(aes(label = value), position = position_stack(), 
  #           vjust = -0.2, hjust = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  labs(x = "", y = "# of large hypo CG-DMRs") +
  theme(legend.title = element_blank(),
        legend.position = "top")
dev.off()

inputs %>% 
  select(c("V1", "V2", "group")) %>% 
  group_by(group) %>% 
  summarise_all(sum) %>%
  filter(group %in% c("Adipose", "Spleen", "Lung", 
                      "Muscle", "Rumen", "Forebrain")) %>% 
  mutate(percent = V1 / V2)
