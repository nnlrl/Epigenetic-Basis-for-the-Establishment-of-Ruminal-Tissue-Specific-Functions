library(gprofiler2)
library(Vennerable)
library(tidyverse)
library(org.Bt.eg.db)
library(clusterProfiler)
library(ChIPseeker)
library(rtracklayer)
library(ggridges)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)


data <- as.data.frame(data.table::fread("./dmr.csv"))
enhancer <- makeGRangesFromDataFrame(
  data.table::fread("./methylation/enhancer/all_enhancer.bed"),
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3"
)

sample.list <- sub("methylation_level_", "", names(data)[8:ncol(data)])

res <- data.frame()
for (sample in c("forebrain", "heart", "muscle", 
                 "lung", "rumen", "kidney")) {
  data.up <- data[grepl(sample, str_to_lower(data$hypermethylated_samples)),
                  c(1:3, grep(sample, str_to_lower(names(data))))]
  data.up$`#chr` <- paste0("chr", data.up$`#chr`)
  up.gr <- makeGRangesFromDataFrame(data.up,
                                    seqnames.field = "#chr",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)
  
  up.2kb <- reduce(up.gr, min.gapwidth = 1000)
  up.2kb <- up.2kb[width(up.2kb) > 2000]
  
  data.down <- data[grepl(sample, str_to_lower(data$hypomethylated_samples)),
                  c(1:3, grep(sample, str_to_lower(names(data))))]
  # data.down$`#chr` <- paste0("chr", data.down$`#chr`)
  down.gr <- makeGRangesFromDataFrame(data.down,
                                    seqnames.field = "#chr",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)
  
  down.2kb <- reduce(down.gr, min.gapwidth = 1000)
  down.2kb <- down.2kb[width(down.2kb) > 2000]
  
  res <- rbind(res, 
               as.data.frame(down.2kb) %>% 
                 mutate(group = sample))
  res$seqnames <- sub("chr", "", res$seqnames)
  data.table::fwrite(as.data.frame(down.2kb)[,c(1,2,3)], paste0(sample, ".2kb.bed"),
                     quote = F, sep = "\t", row.names = F, col.names = F)
}
res$seqnames <- paste0("chr", res$seqnames)
res <- makeGRangesFromDataFrame(res,
                                keep.extra.columns = T)
res$enhancer <- countOverlaps(res, enhancer)

lapply(split(as.data.frame(res), as.data.frame(res)$group),
       function(x){table(x$enhancer)})

pdf("enhancer_intersect2.pdf", 6, 4)
as.data.frame(res) %>% 
  filter(enhancer %in% c(0, 1)) %>% 
  mutate(enhancer = factor(as.character(enhancer))) %>% 
  group_by(group, enhancer) %>% 
  summarise(all = n()) %>% 
  ggplot(aes(x = factor(group, levels = c("kidney", "forebrain", "lung",
                                          "heart", "muscle", "rumen")), 
             y = all, fill=enhancer)) +
  geom_col(position = "stack", width = 0.7) +
  # geom_text(aes(label=value, y=value), vjust = -0.5, position = "stack",
  #           size = 5) +
  # coord_flip() +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("1"="#377EB8", "0"="#E41A1C")) +
  labs(x = "", y = "# of large hypo CG-DMRs") +
  theme(legend.title = element_blank(),
        legend.position = "top")
dev.off()

peakAnno <- as.data.frame(annotatePeak(res, flankDistance = 5000,
                                       TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene))

peaks_up <- list()
peaks_down <- list()
for (sample in c("forebrain", "heart", "muscle", "lung", "rumen", "kidney")) {
  name <- paste0("CGlevel_", sample, ".tsv")
  data <- na.omit(as.data.frame(data.table::fread(name)))
  data_up <- data[data[,ncol(data)]-data[,ncol(data)-1]>0.1,] %>% 
    mutate("Status" = "UP")
  data_down <- data[data[,ncol(data)-1]-data[,ncol(data)]>0.1,] %>% 
    mutate("Status" = "DOWN")
  inputs <- rbind(data_up, data_down)
  inputs$chr <- paste0("chr", inputs$chr)
  
  tmp <- makeGRangesFromDataFrame(inputs, keep.extra.columns = T)
  anno <- annotatePeak(tmp, TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene,
                       flankDistance = 5000, annoDb = "org.Bt.eg.db")
  
  data.table::fwrite(as.data.frame(anno), paste0("annoPeak_", sample, ".tsv"), quote = F,
                     col.names = T, row.names = F, sep = "\t")
  
  # anno_row <- data.frame(
  #   row.names = rownames(inputs),
  #   group = c(rep("up", nrow(data_up)),
  #             rep("down", nrow(data_down))))
  # anno_colors <- list(group=c("up"="#E41A1C", "down"="#377EB8"))
  # # pdf(paste0(sample, ".2kb.heatmap.pdf"), 3, 5)
  # # pheatmap::pheatmap(inputs[,(ncol(inputs)-1):ncol(inputs)],
  # #                    cluster_rows = F, cluster_cols = F,
  # #                    show_rownames = F,
  # #                    cellwidth = 15, annotation_row = anno_row,
  # #                    annotation_colors = anno_colors,
  # #                    labels_col = c(paste0("fetal_", str_to_title(sample)),
  # #                                   paste0("adult_", str_to_title(sample))),
  # #                    angle_col = 45,
  # #                    color = colorRampPalette(c("#334fa2", "white", "#ed2224"))(100))
  # # dev.off()
  # data_up$chr <- paste0("chr", data_up$chr)
  # data_up <- makeGRangesFromDataFrame(data_up, 
  #                                     seqnames.field = "chr",
  #                                     start.field = "start",
  #                                     end.field = "end",
  #                                     keep.extra.columns = T)
  # peaks_up[[sample]] <- data_up
  # 
  # data_down$chr <- paste0("chr", data_down$chr)
  # data_down <- makeGRangesFromDataFrame(data_down,
  #                                       seqnames.field = "chr",
  #                                       start.field = "start",
  #                                       end.field = "end", keep.extra.columns = T)
  # peaks_down[[sample]] <- data_down
  # data$diff <- data[,ncol(data)]-data[,ncol(data)-1]
  # 
  # gene_up <- seq2gene(data_up, tssRegion = c(-3000, 3000), flankDistance = 5000,
  #                     TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene)
  # enrich_up <- gost(gene_up, organism = "btaurus")$result
  # # data.table::fwrite(enrich_up[,1:13], paste0(sample, ".enrich_up.tsv"), quote = F,
  # #                    row.names = F, sep = "\t")
  # gene_down <- seq2gene(data_down, tssRegion = c(-3000, 3000), flankDistance = 5000,
  #                       TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene)
  # enrich_down <- gost(gene_down, organism = "btaurus")$result
  # # data.table::fwrite(enrich_down[,1:13], paste0(sample, ".enrich_down.tsv"), quote = F,
  # #                    row.names = F, sep = "\t")
}

# genes expr
load("./methylation/all_RNA.Rdata")
data[as.data.frame(peakAnnoList$forebrain)$SYMBOL,]
data <- as.data.frame(limma::normalizeBetweenArrays(data))
down <- bitr(seq2gene(peaks_down$forebrain, tssRegion = c(-3000, 3000),
                      flankDistance = 5000, TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene), fromType = "ENTREZID", toType = "SYMBOL",
             OrgDb = org.Bt.eg.db)$SYMBOL
up <- bitr(seq2gene(peaks_up$rumen, tssRegion = c(-3000, 3000),
                    flankDistance = 5000, TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene), 
           fromType = "ENTREZID", toType = "SYMBOL",
           OrgDb = org.Bt.eg.db)$SYMBOL

expr_up <- na.omit(data[as.data.frame(peakAnnoList$forebrain)$SYMBOL, allTraits$group %in% c("forebrain", 
                                                   "fetal_Forebrain")]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(sample = c(rep("forebrain", 4),
                    rep("forebrain", 2)),
         group = c(rep("adult", 4), rep("fetal", 2))) %>% 
  group_by(sample, group) %>% 
  summarise_all(mean) %>% 
  reshape2::melt() %>% 
  mutate(status = "up")

expr_down <- na.omit(data[down, allTraits$group %in% c("Rumen", "fetal_Rumen")]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(sample = c(rep("Rumen", 4),
                    rep("Rumen", 2)),
         group = c(rep("adult", 4), rep("fetal", 2))) %>% 
  group_by(sample, group) %>% 
  summarise_all(mean) %>% 
  reshape2::melt() %>% 
  mutate(status = "down")

rbind(expr_up) %>% 
  ggplot(aes(x = group, y = log10(value + 1) , fill = group)) +
  geom_boxplot()



inputs <- data.frame()
peakAnnoList <- lapply(peaks_down, annotatePeak, 
                       TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene,
                       tssRegion=c(-3000, 3000), flankDistance = 5000)
inputs <- rbind(inputs, 
                unlist(lapply(peakAnnoList, function(i) 
                  length(unique(as.data.frame(i)[as.data.frame(i)$annotation != "Distal Intergenic",]$geneId)))
                ))
pdf("annoBar_down.pdf", 6, 4)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("distToTSS_down.pdf", 6, 4)
plotDistToTSS(peakAnnoList)
dev.off()
genes_down <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)[as.data.frame(i)$annotation != "Distal Intergenic",]$geneId)

peakAnnoList <- lapply(peaks_up, annotatePeak, 
                       TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene,
                       tssRegion=c(-3000, 3000), annoDb = "org.Bt.eg.db")
inputs <- rbind(inputs, 
                unlist(lapply(peakAnnoList, function(i) 
                  length(unique(as.data.frame(i)[as.data.frame(i)$annotation != "Distal Intergenic",]$geneId)))
                ))
names(inputs) <- c("forebrain", "heart", "muscle", "lung", "rumen", "kidney")
inputs$group <- c("down", "up")
pdf("number_of_large_DMRs_Genes.pdf", 4, 3)
inputs %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = factor(variable,
                        levels = c("rumen", "muscle", "heart", "lung", 
                                   "forebrain",
                                   "kidney" 
                                   )), y = value, fill = group)) +
  geom_col(width = 0.6, color="black") +
  scale_fill_manual(values = c('down'="#E41A1C", 'up'='#377EB8'),
                    labels = c("Loss of mCG", 'Gain of mCG')) +
  labs(x="", y="Number of Genes") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  scale_y_reverse()
dev.off()

pdf("annoBar_up.pdf", 6, 4)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("distToTSS_up.pdf", 6, 4)
plotDistToTSS(peakAnnoList)
dev.off()
genes_up <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)[as.data.frame(i)$annotation != "Distal Intergenic",]$geneId)

genes <- mapply(c, genes_up, genes_down, SIMPLIFY=FALSE)

all_genes <- unlist(genes)




library(venn)
venn(genes[c(1,2,3,4,5,6)],
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 0.5,     # 数字大小
     sncs = 1        # 组名字体大小
)
plot(Vennerable::Venn(genes[c(1,2,4)]),  doWeights = TRUE)


inputs.df <- matrix(0, nrow = length(unique(do.call(c, genes))),
                    ncol = 6)
rownames(inputs.df) <- unique(do.call(c, genes))

for (i in 1:length(genes)) {
  for (j in 1:length(genes[[i]])) {
    inputs.df[genes[[i]][j], i] <- 1
  }
}

inputs.df <- as.data.frame(inputs.df)
names(inputs.df) <- names(genes)

library(UpSetR)
pdf("upset.pdf", 6, 5)
upset(inputs.df, nsets = 6, 
      # order.by = "freq",
      mb.ratio = c(0.7, 0.3), nintersects = 20,
      text.scale = 1.5)
dev.off()

tmp <- lapply(peakAnnoList)


genes=lapply(peaks_up, function(i) 
  seq2gene(i, c(-3000, 3000), 5000, TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene))
cc = compareCluster(geneClusters = genes, 
                    fun="enrichGO", ont="ALL", OrgDb = org.Bt.eg.db,
                    pvalueCutoff = 1, qvalueCutoff = 1)
pdf("compareCLuster_up.pdf", 10, 10)
dotplot(cc, showCategory=2)
dev.off()

genes=lapply(peaks_down, function(i) 
  seq2gene(i, c(-3000, 3000), 5000, TxDb=TxDb.Btaurus.UCSC.bosTau9.refGene))
cc = compareCluster(geneClusters = genes, 
                    fun="enrichGO", ont="ALL", OrgDb = org.Bt.eg.db,
                    pvalueCutoff = 1, qvalueCutoff = 1)
pdf("compareCLuster_down.pdf", 10, 10)
dotplot(cc, showCategory=2)
dev.off()

pdf("number_of_large_DMRs.pdf", 4, 3)
rbind(unlist(lapply(peaks_up, length)),
      unlist(lapply(peaks_down, length))) %>% 
  reshape2::melt() %>%
  ggplot(aes(x=factor(str_to_title(Var2), levels = c("Rumen", "Muscle", "Heart",
                                                      "Lung", "Kidney", "Forebrain")),
             y=value,fill=factor(Var1, levels = c('2', '1')))) +
  geom_col(color="black", width = 0.6) +
  scale_fill_manual(values = c('2'="#E41A1C", '1'='#377EB8'),
                    labels = c("Loss of mCG", 'Gain of mCG')) +
  labs(x="", y="Number of large CG-DMRs") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  scale_y_reverse()
dev.off()

# vs WGCNA
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


for (sample in c("forebrain", "rumen", "lung", "heart", "muscle", "kidney")) {
  meth <- as.data.frame(data.table::fread(
    paste0("CGlevel_", sample, ".tsv")))
  meth$chr <- paste0("chr", meth$chr)
  data_up <- meth[meth[,ncol(meth)]-meth[,ncol(meth)-1]>0.1,] %>% 
    mutate(group = "up")
  data_down <- meth[meth[,ncol(meth)-1]-meth[,ncol(meth)]>0.1,] %>% 
    mutate(group = "down")
  inputs <- na.omit(reshape2::melt(rbind(data_up, data_down)[,4:6]))
  
  p <- ggplot(inputs, aes(x=group, y=value, fill=variable)) +
    geom_boxplot(width = 0.5) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1),
                       labels = c(0, 0.25, "0.50", 0.75, 1)) +
    guides(fill=F) +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#344ea2", "#ed2224")) +
    labs(x = "", y = "mCG level", 
         title = 
           paste0(nrow(data_down), " vs ", 
                  nrow(data_up)))
  ggsave(paste0(sample, "_high_vs_low.boxplot.pdf"), p, width = 3, height = 3)
  
}
