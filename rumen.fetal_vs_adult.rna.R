library(tidyverse)
library(edgeR)
library(org.Bt.eg.db)
library(clusterProfiler)
library(gprofiler2)


fetal <- data.table::fread("./methylation/fetal_counts.tsv")[,c(1, 6,7)]
names(fetal) <- c("id", "fetal_Rumen_1", "fetal_Rumen_2")
adult <- data.table::fread("./methylation/adult_counts.tsv")[,c(1, 20, 22, 23, 24)]
names(adult) <- c("id", "adult_Rumen_1", "adult_Rumen_2", "adult_Rumen_3",
                  "adult_Rumen_4")

tmp <- inner_join(fetal, adult, by="id")
data <- inner_join(
  bitr(tmp$id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db),
  tmp,
  by=c("SYMBOL"="id")
)

y <- DGEList(data[,3:ncol(data)], genes = data[,1:2], group = c(rep("fetal", 2),
                                                                rep("adult", 4)))
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)

res <- exactTest(y, pair = c("fetal", "adult"))
diff <- topTags(res, n = Inf)$table

pval <- 0.05
logfc <- 1.5
write.csv(diff, file = "allDiff.csv", quote = FALSE)
diffSig <- diff[(diff$PValue < pval & (diff$logFC > logfc | diff$logFC < -logfc)), ]
write.csv(diffSig, file = "diffSig.csv", quote = FALSE)
diffUp <- diff[(diff$PValue < pval & (diff$logFC > logfc)), ] # for go, kegg
write.csv(diffUp, file = "diffUp.csv", quote = FALSE)
diffDown <- diff[(diff$PValue < pval & (diff$logFC < -logfc)), ] # for go, kegg
write.csv(diffDown, file = "diffDown.csv", quote = FALSE)

res <- gost(diffSig$ENTREZID, organism = "btaurus")$result
data.table::fwrite(res[,1:13], "diffSig.gost.tsv", quote = F,
                   sep = "\t")
res <- gost(diffUp$ENTREZID, organism = "btaurus")$result
data.table::fwrite(res[,1:13], "diffUp.gost.tsv", quote = F,
                   sep = "\t")
res <- gost(diffDown$ENTREZID, organism = "btaurus")$result
data.table::fwrite(res[,1:13], "diffDown.gost.tsv", quote = F,
                   sep = "\t")

go <- as.data.frame(enrichGO(diffDown$SYMBOL, OrgDb = org.Bt.eg.db,
                             keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.1,
                             qvalueCutoff = 1))

data.table::fwrite(go, "go.down.tsv", quote = F, sep = "\t")

anno_go(diffSig$ENTREZID, "all", org.Bt.eg.db, "SYMBOL")
anno_go(diffUp$ENTREZID, "up", org.Bt.eg.db, "SYMBOL")
anno_go(diffDown$ENTREZID, "down", org.Bt.eg.db, "SYMBOL")
anno_kegg(diffSig$ENTREZID, "bovine", "all", org.Bt.eg.db, "ENTREZID")
anno_kegg(diffUp$ENTREZID, "bovine", "up", org.Bt.eg.db, "ENTREZID")
anno_kegg(diffDown$ENTREZID, "bovine", "down", org.Bt.eg.db, "ENTREZID")

inputs <- diff
inputs$SYMBOL[which(!inputs$SYMBOL %in% c(
  "KRT14", "KRT17", "MAPK1", "IGF2BP1", "SOX10", "S100A8", "TGM3", "DLK1", "GAP43", "COL11A1",
  "HMGCS2", "SLC16A1", "NOD2", "DEFB1", "LYZ1", "LYZ2B", "IL17A", "IL17F",
  "IL17RC", "F2RL1", "TMPRSS15"))] <- NA
pdf("volcano.pdf", 10, 10)
p <- EnhancedVolcano::EnhancedVolcano(
  toptable = inputs,
  lab = inputs$SYMBOL,
  x = "logFC", y = "PValue", pCutoff = 0.05, FCcutoff = 1.5,
  subtitle = "fetal vs adult", caption = "",
  ylim = c(0, 200), xlim = c(-20, 20)
)
print(p)
dev.off()

anno_go <- function(y, kind, org, keytype) {
  go_ <- function(data_id, ont, kind, pval, qval) {
    if (!(ont %in% c("BP", "CC", "MF", "ALL"))) stop("Invalid input: ont")
    if (!(str_to_lower(kind) %in% c("up", "down", "all"))) stop("Invalid input: kind")
    data_id <- na.omit(data_id)
    erich.go <- enrichGO(
      gene = data_id,
      OrgDb = org,
      keyType = "ENTREZID",
      ont = ont,
      pvalueCutoff = pval,
      qvalueCutoff = qval
    )
    if (nrow(as.data.frame(erich.go)) > 0) {
      file_name <- paste0("go_", ont, "_", kind)
      write.csv(as.data.frame(erich.go), file = paste0(file_name, ".csv"), quote = F)
      # dotplot
      pdf(file.path(paste0(file_name, "_dotplot.pdf")), 10, 10)
      p <- suppressMessages(dotplot(erich.go, showCategory = 10, font.size = 15))
      print(p)
      dev.off()
      # barplot
      bar_data <- as.data.frame(erich.go)
      if ("ONTOLOGY" %in% colnames(bar_data)) {
        go_list <- split(bar_data, bar_data$ONTOLOGY)
        if ("BP" %in% names(go_list) && "CC" %in% names(go_list) && "MF" %in% names(go_list)) {
          min_row_length <- 10
          res <- lapply(go_list, function(x) {
            x.order <- order(x$Count, decreasing = T)
            x <- x[x.order, ]
            x[1:min_row_length, ]
          })
          res_df <- as.data.frame(do.call(rbind, res))
          res_df <- res_df[!apply(res_df, 1, anyNA), ]
          res_df$Description <- factor(res_df$Description, levels = rev(res_df$Description))
          pdf(file.path(paste0(file_name, "_barplot.pdf")), 10, 10)
          p <- ggplot(res_df, aes(x = Description, y = Count, fill = ONTOLOGY)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_manual(values = c("#6495ED", "#8FBC8F", "#F4A460")) +
            cowplot::theme_cowplot()
          print(p)
          dev.off()
        }
      } else {
        min_row_length <- 20
        res_df <- bar_data[order(bar_data$Count, decreasing = T), ][1:min_row_length, ]
        res_df <- res_df[!apply(res_df, 1, anyNA), ]
        res_df$Description <- factor(res_df$Description, levels = rev(res_df$Description))
        pdf(file.path(paste0(file_name, "_barplot.pdf")), 10, 10)
        p <- ggplot(res_df, aes(x = Description, y = Count)) +
          geom_bar(stat = "identity", fill = "#8FBC8F") +
          coord_flip() +
          cowplot::theme_cowplot()
        print(p)
        dev.off()
      }
    }
  }
  data_id <- y
  # BP
  go_(data_id, "BP", kind, 0.05, 1)
  # CC
  go_(data_id, "CC", kind, 0.05, 1)
  # MF
  go_(data_id, "MF", kind, 0.05, 1)
  # ALL
  go_(data_id, "ALL", kind, 0.05, 1)
}

anno_kegg <- function(y, species, kind, org, keytype) {
  kegg_ <- function(data_id, ont, kind, pval, qval) {
    if (!(str_to_lower(kind) %in% c("up", "down", "all"))) stop("Invalid input: kind")
    data_id <- na.omit(data_id)
    erich.kegg <- suppressMessages(enrichKEGG(
      data_id,
      organism = ont,
      keyType = "kegg",
      pvalueCutoff = pval,
      qvalueCutoff = qval,
      pAdjustMethod = "BH"
    ))
    if (nrow(as.data.frame(erich.kegg)) > 0) {
      file_name <- paste0("kegg_", kind)
      write.csv(as.data.frame(erich.kegg), file = paste0(file_name, ".csv"), quote = F)
      # dotplot
      pdf(file.path(paste0(file_name, "_dotplot.pdf")), 10, 10)
      p <- suppressMessages(dotplot(erich.kegg, showCategory = 10, font.size = 15))
      print(p)
      dev.off()
      
      # barplot
      bar_data <- as.data.frame(erich.kegg)
      min_row_length <- 20
      res_df <- bar_data[order(bar_data$Count, decreasing = T), ][1:min_row_length, ]
      res_df <- res_df[!apply(res_df, 1, anyNA), ]
      res_df$Description <- factor(res_df$Description, levels = rev(res_df$Description))
      pdf(file.path(paste0(file_name, "_barplot.pdf")), 10, 10)
      p <- ggplot(res_df, aes(x = Description, y = Count)) +
        geom_bar(stat = "identity", fill = "#8FBC8F") +
        coord_flip() +
        cowplot::theme_cowplot()
      print(p)
      dev.off()
    }
  }
  
  ont_list <- list(
    human = "hsa",
    bovine = "bta",
    mouse = "mmu",
    goat = "chx",
    sheep = "oas"
  )
  data_id <- y
  kegg_(data_id, ont_list[[species]], kind, 0.05, 1)
}
