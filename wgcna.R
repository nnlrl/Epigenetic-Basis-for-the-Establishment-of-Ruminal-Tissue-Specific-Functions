library(WGCNA)
library(tidyverse)
library(dendextend)
library(org.Bt.eg.db)
library(clusterProfiler)


stemCells <- data.table::fread("./methylation/stemCells.tpm.tsv")
tissues <- data.table::fread("./methylation/GSE137943_UMD3.1.1.fpkm.txt")
genes <- bitr(stemCells$V1, fromType = "SYMBOL", toType = "ENSEMBL",
              OrgDb = org.Bt.eg.db)
genes <- genes[!duplicated(genes$SYMBOL),]

data <- inner_join(inner_join(genes, stemCells, by=c("SYMBOL"="V1")),
                   tissues, by=c("ENSEMBL"="V1")) %>% 
  remove_rownames() %>% 
  dplyr::select(-c("ENSEMBL")) %>% 
  column_to_rownames(var = "SYMBOL")

rawData <- data
data <- log2(data+1) %>% 
  t() %>% 
  as.data.frame()

allTraits <- data.frame(
  row.names = rownames(data),
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

out <- !(allTraits$group %in% c("Lymphocyte", "Ovary"))
data <- data[out,]
allTraits <- allTraits %>% 
  filter(out)

identical(rownames(data), rownames(allTraits))
my_mad <- function(x){mad(x,na.rm = TRUE)}
m.mad <- apply(data,2,my_mad)

# datExpr <- datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:2500],]
datExpr <- data[,which(m.mad > quantile(m.mad, probs=seq(0, 1, 0.25))[2])]
rawData <- rawData[match(colnames(datExpr), rownames(rawData)),]
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = as.dendrogram(hclust(dist(datExpr), method = "ward.D2"))
my_colors <- ifelse(startsWith(allTraits$group, 'b'), "forestgreen", "green")
# sizeGrWindow(1,1) #视图
pdf("sampleCluster.pdf", 12, 5)
par(cex = 0.8);
par(mar = c(8, 1, 2, 1) + 0.1)
sampleTree %>% 
  set("labels_col", value = c("skyblue", "orange"), k=2) %>% 
  set("branches_k_color", value = c("skyblue", "orange"), k = 2) %>%
  plot(axes=FALSE)
# rect.dendrogram(sampleTree, k=3, lty = 5, lwd = 0, x=1, col=rgb(0.1, 0.2, 0.4, 0.1) )
# colored_bars(colors = my_colors, dend = sampleTree, rowLabels = "am")
dev.off()


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

if (T) {
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.9,col="red")                  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
sft$powerEstimate
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)
cor <- stats::cor

net <- readRDS("./methylation/wgcna/net.Rds")

table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
p <- pheatmap::pheatmap(MEs, clustering_method = "ward.D2",
                        scale = "row",
                        cellwidth = 10,
                        border_color = NA,
                        colorRampPalette(c("#334fa2", "white", "#ed2224"))(100))

MEs <- MEs[p$tree_row$order,p$tree_col$order]
MEs <- MEs[c(1:34, 36:nrow(MEs), 35),]
pdf("allMEs_heatmap.pdf", 6, 8)
pheatmap::pheatmap(MEs, clustering_method = "ward.D2",
                   scale = "row",
                   cluster_rows = F,
                   cluster_cols = F,
                   cellwidth = 10,
                   border_color = NA,
                   color = colorRampPalette(c("#334fa2", "white", "#ed2224"))(100),
                   gaps_row = c(3, 18, 22, 26, 30, 34, 36, 40, 44, 48, 52,
                                56, 60, 62))
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

for (module in unique(moduleColors)) {
  tmp <- moduleColors == module
  inputs <- rawData[tmp,]
  data.table::fwrite(inputs, paste0(module, "_sub.tpm.tsv"),
                     quote = F, sep = "\t", row.names = T)
  
  # inputs <- as.matrix(log10(inputs+1))
  # pdf(paste0(module, "_sub_heatmap.pdf"), 6, 8)
  # pheatmap::pheatmap(inputs, show_rownames = F, scale = "row",
  #                    cluster_cols = F,
  #                    cluster_rows = F,
  #                    border_color = NA,
  #                    color = colorRampPalette(c("#334fa2", "white", "#ed2224"))(100))
  # dev.off()
  # 
  # inputs <- reshape2::melt(inputs)
  # inputs$Var2 <- sub("_[12345]$", "", inputs$Var2)
  # inputs$Var2 <- sub("bESCs_F7.+", "bESCs_F7", inputs$Var2)
  # inputs$value <- as.numeric(inputs$value)
  # 
  # pdf(paste0(module, "_sub_boxplot.pdf"), 7, 5)
  # p <- inputs %>%
  #   group_by(Var1, Var2) %>%
  #   summarise(x = mean(value)) %>%
  #   ggplot(aes(x=Var2, y=x)) +
  #   geom_boxplot(fill="grey") +
  #   cowplot::theme_cowplot() +
  #   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  #   labs(x="Samples", y=expression(log[10](TPM+1)))
  # print(p)
  # dev.off()
}


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
# Select module

for (module in unique(moduleColors)) {
  # Select module probes
  probes = colnames(datExpr)
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.2,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  )
}

for (module in c("black", "blue", "brown", "cyan", "green",
                 "greenyellow", "grey", "grey60", "lightcyan",
                 "lightgreen", "magenta", "midnightblue",
                 "pink", "purple", "red", "salmon", "tan",
                 "turquoise", "yellow")) {
  file <- paste0(module, "_sub.tpm.tsv")
  
  data <- data.table::fread(paste0("./methylation/wgcna/", file))
  genes <- data$V1
  enrich <- gost(genes, organism = "btaurus")$result
  data.table::fwrite(enrich[,1:13], paste0("./methylation/wgcna/enrich_",
                                           module, ".tsv"),
                     quote = F, sep = "\t", row.names = F, col.names = T)
}

data <- data.table::fread(paste0("./methylation/wgcna/tan_sub.tpm.tsv"))
genes <- data$V1
enrich <- gost(genes, organism = "btaurus")$result

tfs <- data.table::fread("./methylation/wgcna/TF_names_v_1.01.txt", header = F)
intersect(tfs$V1, genes)

res <- enrichGO(data$V1, OrgDb = org.Bt.eg.db, keyType = "SYMBOL",
                ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1)
View(as.data.frame(res))

moduleColors[which(names(datExpr) == "DNMT1")]
