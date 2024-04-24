library(ggplot2)

stitched_regions <- read.delim(
  file= "./methylation/rose/rose/rumen_1/adult_rumen_H3K27ac_REP1_peaks_12.5KB_STITCHED_TSS_DISTAL_REGION_MAP.txt",sep="\t")
rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])
rankBy_vector[rankBy_vector <= 0] <- 0
signalOrder = order(rankBy_vector,decreasing=TRUE)

cutoff <- 57852.9336 # cutoff of rose
stitched_regions$group <- "x"
stitched_regions$group[which(rankBy_vector >= cutoff)] <- "y"
superEnhancerRows <- which(rankBy_vector > cutoff)

# plot(length(rankBy_vector):1,rankBy_vector[signalOrder], 
#      col='red',xlab=paste(rankBy_factor,' Stitched peaks'),
#      ylab=paste(rankBy_factor,' Signal'),pch=19,cex=2)
# abline(h=cutoff,col='grey',lty=2)

pdf("./methylation/rose/fetal_H3K27me3_points2.pdf", 3, 3)
ggplot2::ggplot(data = data.frame()) +
  geom_line(aes(length(rankBy_vector):1, rankBy_vector[signalOrder]),
            color = "grey60", linewidth = 1) +
  geom_point(aes(x = length(rankBy_vector):1,
                 y = rankBy_vector[signalOrder],
                 color = stitched_regions$group[signalOrder]),
             size = 2) +
  geom_hline(aes(yintercept = cutoff),
             color = "grey", linetype = "dashed") +
  geom_vline(aes(xintercept = length(rankBy_vector)-length(superEnhancerRows)),
             color = "grey", linetype = "dashed") +
  cowplot::theme_cowplot() +
  theme(legend.position = 'none') +
  labs(x = "Rank", y = "H3K27me3 signal")
dev.off()


# convert gene name
tmp <- read.delim(
  file= "./methylation/rose/rose/adult_H3K27me3_rumen_1/adult_rumen_H3K27me3_REP1_peaks_AllStitched_REGION_TO_GENE.txt",sep="\t") %>% 
  filter(isSuper == 1)

genes <- bitr(
  unique(c(
    tmp$CLOSEST_GENE, 
    unlist(strsplit(tmp$OVERLAP_GENES, split = ",")),
    unlist(strsplit(tmp$PROXIMAL_GENES, split = ","))
  )), fromType = "REFSEQ", toType = "SYMBOL",
  OrgDb = org.Bt.eg.db
)

tmp <- t(apply(tmp, 1, function(x) {
  for (i in 1:nrow(genes)) {
    x <- gsub(genes[i, 1], genes[i, 2], x)
  }
  x
}))

data.table::fwrite(tmp, 
                   "./methylation/rose/rose/adult_H3K27me3_rumen_1/adult_rumen_H3K27me3_REP1_peaks_AllStitched_REGION_TO_GENE.SYMBOL.txt",
                   quote = F, sep = "\t", row.names = F, col.names = T)

library(tidyverse)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
enhancer <- GenomicRanges::union(
  makeGRangesFromDataFrame(
    read.delim(
      file= "./methylation/rose/rose/adult_H3K27me3_rumen_1/adult_rumen_H3K27me3_REP1_peaks_AllStitched_REGION_TO_GENE.txt",sep="\t") %>% 
      filter(isSuper == 1) %>% 
      select(c("CHROM", "START", "STOP")),
    seqnames.field = "CHROM",
    start.field = "START",
    end.field = "STOP"
  ), makeGRangesFromDataFrame(
    read.delim(
      file= "./methylation/rose/rose/adult_H3K27me3_rumen_2/adult_rumen_H3K27me3_REP2_peaks_AllStitched_REGION_TO_GENE.txt",sep="\t") %>% 
      filter(isSuper == 1) %>% 
      select(c("CHROM", "START", "STOP")),
    seqnames.field = "CHROM",
    start.field = "START",
    end.field = "STOP"
  )
)


dmrs <- data.table::fread("D:/METH/rumen/dmrs.tsv")
dmrs$chr <- paste0("chr", dmrs$chr)
dmrs <- makeGRangesFromDataFrame(dmrs, keep.extra.columns = T)

feDMRs <- as.data.frame(dmrs[countOverlaps(dmrs, enhancer) > 0])

data.table::fwrite(feDMRs, file = "adult_rumen.H3K27me3.bed", quote = F,
                   sep = "\t", row.names = F, col.names = F)


# barplot
inputs <- data.frame(
  group = c("fetal", "adult"),
  nums = c(2557, 92)
)

pdf("number_silencer.pdf", 3, 5)
ggplot(inputs, aes(x = factor(group, levels = c("fetal", "adult")),
                   y = nums, fill = group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Stage", y = "Number") +
  cowplot::theme_cowplot()
dev.off()


# chip & rna
mCG <- data.table::fread("./methylation/rose/RefSeqPromoter.Rumen.tsv")
allDiff <- data.table::fread("./methylation/rose/allDiff.csv")
adult_H3K27ac <- data.table::fread("./methylation/rose/rose/rumen_2/adult_rumen_H3K27ac_REP2_peaks_SuperStitched_GENE_TO_REGION.txt")
fetal_H3K27ac <- data.table::fread("./methylation/rose/rose/fetal_rumen_4/fetal_rumen_H3K27ac_REP2_peaks_SuperStitched_GENE_TO_REGION.txt")
adult_H3K27me3 <- data.table::fread("./methylation/rose/rose/adult_H3K27me3_rumen_1/adult_rumen_H3K27me3_REP1_peaks_SuperStitched_GENE_TO_REGION.txt")
fetal_H3K27me3 <- data.table::fread("./methylation/rose/rose/fetal_H3K27me3_rumen_2/fetal_rumen_H3K27me3_REP2_peaks_SuperStitched_GENE_TO_REGION.txt")


data <- left_join(
  left_join(
    left_join(
      left_join(
        inner_join(
          mCG, allDiff, by = c("col_3" = "ENTREZID")
        ), adult_H3K27ac,
        by = c("SYMBOL" = "GENE_NAME")
      ), fetal_H3K27ac, by = c("SYMBOL" = "GENE_NAME")
    ), adult_H3K27me3, by = c("SYMBOL" = "GENE_NAME")
  ), fetal_H3K27me3, by = c("SYMBOL" = "GENE_NAME")
)

data$adult_H3K27ac <- ifelse(is.na(data$PROXIMAL_STITCHED_PEAKS.x), "no", "yes")
data$fetal_H3K27ac <- ifelse(is.na(data$PROXIMAL_STITCHED_PEAKS.y), "no", "yes")
data$adult_H3K27me3 <- ifelse(is.na(data$PROXIMAL_STITCHED_PEAKS.x.x), "no", "yes")
data$fetal_H3K27me3 <- ifelse(is.na(data$PROXIMAL_STITCHED_PEAKS.y.y), "no", "yes")
data$methylation_level_adult <- apply(data[,c(6,7)], 1, mean)
data$methDiff <- data$methylation_level_adult - data$methylation_level_fetal
data$logFC[data$logFC > 10] <- 10
tmp <- data[data$logFC > 2 & data$methDiff < -0.2 & data$adult_H3K27ac == "yes"]
tmp2 <- data[data$adult_H3K27me3 == "yes"]

data$repel <- data$SYMBOL
data$repel[data$adult_H3K27me3 == "no" & (data$logFC < 2 | data$methDiff > -0.2 | data$adult_H3K27ac == "no")] <- NA

data$part <- case_when(abs(data$logFC) >= 1.5 & abs(data$methDiff) >= 0.2 ~ "part1379",
                       data$logFC >= 1.5 & data$methDiff >= 0.2 ~ "part3",
                       abs(data$logFC) < 1.5 & abs(data$methDiff) > 0.2 ~ "part28",
                       abs(data$logFC) > 1.5 & abs(data$methDiff) < 0.2 ~ "part46",
                       abs(data$logFC) < 1.5 & abs(data$methDiff) < 0.2 ~ "part5")
head(data)

pdf("RNA&mCG.points.2.pdf", 7, 7)
ggplot(data, aes(x = logFC, 
                 y = methylation_level_adult - methylation_level_fetal,
                 color = part)) +
  geom_point() +
  geom_text_repel(aes(label = repel), color = "black") +
  cowplot::theme_cowplot() +
  ylim(c(-1, 1)) +
  xlim(c(-10, 10)) +
  scale_colour_manual(name="",values=alpha(c("#872d25", "#eeea33", "#6e76b2",  "gray80"),0.7)) +
  geom_hline(yintercept = c(-0.2,0.2),
             size = 0.5,
             color = "grey40",
             lty = "dashed")+
  geom_vline(xintercept = c(-1.5,1.5),
             size = 0.5,
             color = "grey40",
             lty = "dashed") +
  geom_point(data = tmp, aes(x = logFC, y = methDiff), color = "#4997c7") +
  geom_point(data = tmp2, aes(x = logFC, y = methDiff), color = "#8634a6")
  
dev.off()


data$group <- "other"
data$group[data$methDiff >= 0.2 | data$methDiff <= -0.2] <- "mCG"
data$group[(data$adult_H3K27ac == "yes" | data$fetal_H3K27ac == "yes") &
             data$fetal_H3K27me3 == "no" & data$adult_H3K27me3 == "no"] <- "H3K27ac"
data$group[data$adult_H3K27ac == "no" & data$fetal_H3K27ac == "no" &
             (data$fetal_H3K27me3 == "yes" | data$adult_H3K27me3 == "yes")] <- "H3K27me3"

diffUp <- data[data$logFC >= 1.5,]
diffDown <- data[data$logFC <= -1.5,]
pie(table(diffUp$group))
pie(table(diffDown$group))

# gsea
library(clusterProfiler)

H3K27ac <- data[data$adult_H3K27ac == "yes" ,][,c(4,10)]
H3K27ac <- H3K27ac[order(H3K27ac$logFC, decreasing = T),]

H3K27me3 <- data[data$fetal_H3K27ac == "yes",][,c(4,10)]
H3K27me3 <- H3K27me3[order(H3K27me3$logFC, decreasing = T),]

res <- enrichKEGG(H3K27ac$col_3, organism = 'bta', minGSSize = 3, maxGSSize = 1000, pvalueCutoff=1)
View(as.data.frame(res))
pdf("adult_H3K27ac.pdf", 7, 5)
enrichplot::gseaplot2(res, 1:5, pvalue_table = TRUE)
dev.off()

res_df <- as.data.frame(res)[c(3,4,5,6,9),]
res_df <- res_df[order(res_df$Count, decreasing = T), ]
res_df$Description <- factor(res_df$Description, levels = rev(res_df$Description))
pdf(paste0("fetal_H3K27ac_barplot.pdf"), 10, 10)
p <- ggplot(res_df, aes(x = Description, y = Count)) +
  geom_bar(stat = "identity", fill = "grey60") +
  coord_flip() +
  cowplot::theme_cowplot()
print(p)
dev.off()
