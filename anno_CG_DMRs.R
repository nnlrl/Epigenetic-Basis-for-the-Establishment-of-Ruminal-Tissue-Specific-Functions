library(cowplot)
library(tidyverse)
library(ChIPseeker)
library(org.Bt.eg.db)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(rtracklayer)
library(clusterProfiler)
library(GenomicFeatures)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
setwd("D:/METH/")

txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene

data <- data.table::fread("./dmr.csv")
# data <- data[,2:ncol(data)]
data$`#chr` <- paste0("chr", data$`#chr`)
data <- data[data$hypermethylated_samples != "" | 
               data$hypomethylated_samples != "",]

# tmp <- apply(data[,8:ncol(data)], 1, function(x) {
#  mean(na.omit(x))
# })

data <- makeGRangesFromDataFrame(
  data,
  keep.extra.columns = TRUE,
  seqnames.field = "#chr",
  start.field = "start",
  end.field = "end"
)

data$length <- width(data)
mean_len <- mean(data$length)
# data$length[data$length > 2000] <- 2001
  
ggplot(as.data.frame(data), aes(length)) +
  geom_histogram(binwidth = 200, color = "white", size=1, alpha=0.7) +
  theme_cowplot() +
  labs(x = "CG-DMR size(bp)", y = "Number of CG-DMRs")

# covplot(data)

annotation.isExistingFileOrNone <- function(fileName) {
  if (is.null(fileName)) {
    return(FALSE)
  }
  
  if (length(fileName) == 1 && file.exists(fileName)) {
    return(TRUE)
  }
  
  stop(paste0("Annotation file not found: ", fileName))
}
annotation.annotateCGIslands <- function(regions, regionRanges, cgis, expend = 2000) {
  cgiRanges <- GRanges(seqnames = cgis$chrom, 
                       ranges = IRanges(cgis$chromStart, cgis$chromEnd))
  cgiShore <- c(flank(cgiRanges, width = expend), 
                flank(cgiRanges, width = expend, start = F))
  regions$cgi_overlap <- countOverlaps(regionRanges, cgiRanges) > 0
  regions$cgi_shore_overlap <- countOverlaps(regionRanges, cgiShore) > 0
  return(regions)
}
annotation.annotateRepeats <- function(regionTable, regionRanges, repeats, classes) {
  for (class in classes) {
    relevantRepeats <- repeats[repClass %in% class]
    repeatRanges <- GRanges(relevantRepeats$genoName, 
                            IRanges(relevantRepeats$genoStart, 
                                    relevantRepeats$genoEnd))
    regionTable[,class] <- countOverlaps(regionRanges, repeatRanges) > 0
  }
  return(regionTable)
}
annotation.annotateEnhancer <- function(regions, regionRanges, 
                                        activeEnhancer, weakEnhancer, respress) {

  regions$active_enhancer <- countOverlaps(regionRanges, activeEnhancer) > 0
  regions$weak_enhancer <- countOverlaps(regionRanges, weakEnhancer) > 0
  regions$respress <- countOverlaps(regionRanges, respress) > 0
  return(regions)
}
annotation.annotateRegions <- function(regionTable, 
                                       cgiFile, 
                                       txdb, 
                                       repeatMaskerAnnotationFile, 
                                       weakEnhancer,
                                       activeEnhancer,
                                       respress,
                                       promoterTSSDistances) {
  
  promoterTSSDistances <- as.numeric(promoterTSSDistances)
  regionTable$length <- regionTable$end - regionTable$start
  regionRanges <- GRanges(seqnames = regionTable$seqnames, 
                          ranges = IRanges(regionTable$start, 
                                           regionTable$end))
  
  if (annotation.isExistingFileOrNone(cgiFile)) {
    message("Annotation CgiFile...")
    cgis <- fread(cgiFile)
    regionTable <- annotation.annotateCGIslands(regionTable, regionRanges, cgis)
  }
  
  if (T) {
    message("Annotation PromoterFile...")
    res <- ChIPseeker::annotatePeak(regionRanges,
                                    tssRegion = c(-promoterTSSDistances, 
                                                  promoterTSSDistances),
                                    TxDb = txdb,
                                    addFlankGeneInfo = TRUE,
                                    flankDistance = 5000
    )
    res <- as.data.frame(res)
    if (identical(res$start, regionTable$start) &&
        identical(res$end, regionTable$end) &&
        identical(res$seqnames, regionTable$seqnames)) {
      regionTable$annotation <- res$annotation
      regionTable$geneId <- res$geneId
      regionTable$distanceToTSS <- res$distanceToTSS
    }
    
    # regionTable <- annotation.annotateGenes(regionTable, regionRanges, gtfRanges)
    # regionTable <- annotation.annotatePromoters(regionTable, regionRanges, gtfRanges, promoterTSSDistances)
  }
  
  
  message("Annotation PromoterFile...")
  regionTable <- annotation.annotateEnhancer(regionTable, regionRanges, 
                                             activeEnhancer, weakEnhancer, respress)
  
  if (annotation.isExistingFileOrNone(repeatMaskerAnnotationFile)) {
    message("Annotation repeatFile...")
    repeats <- fread(repeatMaskerAnnotationFile)
    regionTable <- annotation.annotateRepeats(
      regionTable, 
      regionRanges, 
      repeats, 
      c("LINE", "SINE", "LTR", "DNA", "Simple_repeat", "Satellite"))
  }
  
  return(regionTable)
}

weakEnhancer <- GRanges()
activeEnhancer <- GRanges()
respress <- GRanges()
for (sample in c("Cortex", "Lung", "Muscle", "adult_Rumen", "fetal_Rumen", "Spleen", "Adipose")) {
  file <- paste0("./methylation/chromhmm/learnmodel_10/learnmodel_10/", sample, "_10_chrhmm_segments.bed.gz")
  df <- data.table::fread(file)
  gr <- makeGRangesFromDataFrame(
    df,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    keep.extra.columns = T
  )
  weakEnhancer <- union(weakEnhancer, gr[gr$V4 %in% c("E9", "E10")])  
  activeEnhancer <- union(activeEnhancer, gr[gr$V4 %in% c("E8")])  
  respress <- union(respress, gr[gr$V4 %in% c("E3")])
}



regionTable <- as.data.frame(data)[,1:3]
regionTable <- annotation.annotateRegions(regionTable,
                           "D:/METH/cpg.bed",
                           txdb,
                           "D:/METH/repeat.bed",
                           weakEnhancer = weakEnhancer,
                           activeEnhancer =  activeEnhancer,
                           respress =  respress,
                           promoterTSSDistances = 3000)

regionTable$res <- ""
regionTable$annotation <- trimws(str_split(regionTable$annotation, '[(]', simplify = T)[,1], which = "both")
regionTable$res[regionTable$annotation == "Promoter" & (regionTable$cgi_overlap | regionTable$cgi_shore_overlap)] <- "CGI Promoter"
regionTable$res[regionTable$res == "" & regionTable$annotation == "Promoter"] <- "non-CGI Promoter"
regionTable$res[regionTable$res == "" & regionTable$cgi_overlap] <- "CGI"
regionTable$res[regionTable$res == "" & 
                  regionTable$annotation %in% c("3' UTR", "5' UTR", "Downstream",
                                                "Exon", "Intron")] <- "Intragenic"
regionTable$res[regionTable$res == "" & regionTable$active_enhancer] <- "Active Enhancer"
regionTable$res[regionTable$res == "" & regionTable$weak_enhancer] <- "Weak/Poised Enhancer"
regionTable$res[regionTable$res == "" & regionTable$cgi_shore_overlap] <- "CGI shore"
# regionTable$res[regionTable$res == "" & regionTable$respress] <- "Polycomb-repressed"
regionTable$res[regionTable$res == "" & regionTable$LTR] <- "LTR"
regionTable$res[regionTable$res == "" & regionTable$LINE] <- "LINE"
regionTable$res[regionTable$res == "" & regionTable$SINE] <- "SINE"
regionTable$res[regionTable$res == "" & regionTable$DNA] <- "DNA repeat"
regionTable$res[regionTable$res == "" & regionTable$Satellite] <- "Satellite"
regionTable$res[regionTable$res == "" & regionTable$Simple_repeat] <- "Simple_repeat"
regionTable$res[regionTable$res == ""] <- "Intergenic"

table(regionTable$res)
pie(table(regionTable$res))

tmp <- regionTable[regionTable$res %in% c("CGI", "CGI Promoter", "CGI shore",
                                          "non-CGI Promoter"),]
pie(table(tmp$res))

sanky_list <- list(
  nodes = data.frame(name = c("All", "Proximal", "Distal", 
                              names(table(regionTable$res)),
                              "Repeats", "Enhancer")),
  links = data.frame(
    source = c(0, 0, 18,1, 1, 1, 17, 2, 0, 17, 17, 1,  17, 17, 17, 18, 2, 2),
    target = c(1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    value = c(59871, 683765, 68284, 12476, 20759, 11417, 11024, 141921,
              201761, 153940, 49563, 15219, 1007, 6112, 90153, 166236,
              311799, 234520)
  )
)
library(networkD3)
sankeyNetwork(Links = sanky_list$links, Nodes = sanky_list$nodes, 
              Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 0, nodeWidth = 30)

regionTable$res <- ""
regionTable$annotation <- trimws(str_split(regionTable$annotation, '[(]', simplify = T)[,1], which = "both")
regionTable$res[regionTable$annotation == "Promoter" & (regionTable$cgi_overlap | regionTable$cgi_shore_overlap)] <- 
  "Proximal CG-DMRs\n(6.78%, n = 64346)"
regionTable$res[regionTable$res == "" & regionTable$annotation == "Promoter"] <- "Proximal CG-DMRs\n(6.78%, n = 64346)"
regionTable$res[regionTable$res == "" & regionTable$cgi_overlap] <- "Proximal CG-DMRs\n(6.78%, n = 64346)"
regionTable$res[regionTable$res == "" & 
                  regionTable$annotation %in% c("3' UTR", "5' UTR", "Downstream",
                                                "Exon", "Intron")] <- "Intragenic\n(21.24%, n = 201761)"

regionTable$res[regionTable$res == "" & regionTable$active_enhancer] <- "Active Enhancer\n(6.94%, n = 65874)"
regionTable$res[regionTable$res == "" & regionTable$weak_enhancer] <- "Weak Enhancer\n(11.78%, n = 111857)"
regionTable$res[regionTable$res == "" & regionTable$cgi_shore_overlap] <- "Proximal CG-DMRs\n(6.78%, n = 64346)"

regionTable$res[regionTable$res == "" & regionTable$respress] <- "Polycomb-repressed\n(11.52%, n = 109398)"

regionTable$res[regionTable$res == "" & regionTable$LTR] <- "Repeats\n(28.38%, n = 269601)"
regionTable$res[regionTable$res == "" & regionTable$LINE] <- "Repeats\n(28.38%, n = 269601)"
regionTable$res[regionTable$res == "" & regionTable$SINE] <- "Repeats\n(28.38%, n = 269601)"
regionTable$res[regionTable$res == "" & regionTable$DNA] <- "Repeats\n(28.38%, n = 269601)"
# regionTable$res[regionTable$res == "" & regionTable$Satellite] <- "Satellite"
# regionTable$res[regionTable$res == "" & regionTable$Simple_repeat] <- "Repeats\n(28.38%, n = 269601)"
regionTable$res[regionTable$res == ""] <- "Intergenic\n(13.37%, n = 127035)"
pie(table(regionTable$res))

tmp <- regionTable[regionTable$res == "Proximal CG-DMRs\n(6.78%, n = 64346)",]
tmp$res <- ""
tmp$res[tmp$annotation == "Promoter" & (tmp$cgi_overlap | tmp$cgi_shore_overlap)] <- "CGI Promoter\n(32.26%, n = 20759)"
tmp$res[tmp$res == "" & tmp$annotation == "Promoter"] <- "non-CGI Promoter\n(23.65%, n = 15219)"
tmp$res[tmp$res == "" & tmp$cgi_overlap] <- "CGIs\n(19.39%, n = 12476)"
tmp$res[tmp$res == "" & tmp$cgi_shore_overlap] <- "CGI shore\n(24.70%, n = 15892)"
pie(table(tmp$res))

tmp <- regionTable[regionTable$res == "Repeats\n(28.38%, n = 269601)",]
tmp$res <- ""
tmp$res[tmp$res == "" & tmp$LTR] <- "LTR (15.12%, n = 40753)"
tmp$res[tmp$res == "" & tmp$LINE] <- "LINE (52.55%, n = 141676)"
tmp$res[tmp$res == "" & tmp$SINE] <- "SINE (28.74%, n = 77491)"
tmp$res[tmp$res == "" & tmp$DNA] <- "DNA Repeat (3.60%, n = 9681)"
# regionTable$res[regionTable$res == "" & regionTable$Satellite] <- "Satellite"
# regionTable$res[regionTable$res == "" & regionTable$Simple_repeat] <- "Simple Repeat (1.68%, n = 11188)"
pie(table(tmp$res))

tmp <- reduce(data[grepl("brain", data$hypomethylated_samples)],
              min.gapwidth = 1000)
tmp <- tmp[width(tmp) >= 2000]

brain <- na.omit(data.df[c(grepl("brain", data$hypomethylated_samples), 
                   grepl("brain", data$hypermethylated_samples)),
                 c("methylation_level_adult_forebrain",
                   "methylation_level_fetal_hindbrain",
                   "methylation_level_fetal_forebrain")])
brain <- brain[brain$methylation_level_adult_forebrain - brain$methylation_level_fetal_hindbrain < -0.5,]
brain <- brain[brain$methylation_level_adult_forebrain - brain$methylation_level_fetal_hindbrain > 0.5,]
pheatmap::pheatmap(brain)

tmp.df <- as.data.frame(tmp)
ggplot(tmp.df, aes(x=res, y=length, fill=res)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0,1500)) +
  theme_cowplot() +
  scale_y_continuous(breaks=c(0, 500, 1000, 1500), 
                     labels = c(0, 500, 1000, ">1500")) +
  scale_x_discrete(labels=c("CGI Promoter",
                            "CGI Shore",
                            "CGI",
                            "non CGI Promoter")) +
  labs(x = "Proximal CG-DMRs", y = "CG-DMRs size(bp)")
