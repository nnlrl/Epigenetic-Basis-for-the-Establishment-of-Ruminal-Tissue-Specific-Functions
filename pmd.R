library(tidyverse)
library(ChIPseeker)
library(org.Bt.eg.db)
library(clusterProfiler)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
txdb <- TxDb.Btaurus.UCSC.bosTau9.refGene

peak <- data.table::fread("./methylation/pmd/pmd_res.final.bed")
peak$V1 <- paste0("chr", peak$V1)

peak <- GRanges(
  peak$V1,
  IRanges(peak$V2, peak$V3)
)

x = annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb,
                 annoDb = "org.Bt.eg.db")
View(as.data.frame(x))
plotAnnoPie(x)

genes <- seq2gene(peak, c(-3000, 3000), 5000, TxDb = txdb)

go <- enrichGO(genes, OrgDb = org.Bt.eg.db, ont = "ALL")
dotplot(go)
View(as.data.frame(go))

kegg <- enrichKEGG(genes, organism = "bta")
dotplot(kegg)

View(as.data.frame(kegg))

df_na <- data.frame(
  value = seq(1, 20),
  x = runif(20),
  y = runif(20),
  z1 = c(rep(NA, 10), rnorm(10))
)

ggplot(df_na, aes(x = value, y)) +
  geom_bar(aes(fill = z1), stat = "identity") +
  scale_fill_gradient(low='#6699CC',high='#CC3333') +
  labs(title = "scale_fill_gradient(low='#6699CC',high='#CC3333')")


data <- data.table::fread("./methylation/GSE137943_UMD3.1.1.fpkm.txt.gz") %>% 
  select(c("V1", "E1_LIVE3770", "E1_LIVE3773", "E1_LIVE3842",  
           "E1_LIVE3886"))

fetal <- data.table::fread("./methylation/pmd/fetal_liver.tsv")
gene_length <- data.table::fread("./methylation/pmd/fetal_liver_2.read_counts.csv")

fetal <- inner_join(gene_length[,c(1,6)], fetal,
                    by=c("Geneid"="V1")) %>% 
  remove_rownames() %>% 
  column_to_rownames("Geneid")



fetal <- fetal %>% 
  rownames_to_column("SYMBOL")

tmp <- bitr(fetal$SYMBOL, fromType = "SYMBOL", toType = "ENSEMBL",
            OrgDb = org.Bt.eg.db)

res <- inner_join(tmp, fetal, by="SYMBOL")

data <- inner_join(res, data, by=c("ENSEMBL"="V1"))


genes <- bitr(seq2gene(peak, c(-3000, 3000), 0, TxDb = txdb),
              fromType = "ENTREZID", toType = "ENSEMBL", 
              OrgDb = org.Bt.eg.db)$ENSEMBL

up_down_100kb <- c(flank(peak, 100000, start = F),
                   flank(peak, 100000, start = T))

genes2 <- bitr(seq2gene(up_down_100kb, c(-3000, 3000), 0, TxDb = txdb),
               fromType = "ENTREZID", toType = "ENSEMBL", 
               OrgDb = org.Bt.eg.db)$ENSEMBL
genes3 <- genes2[!genes2 %in% genes]

tmp <- bitr(genes3,
            fromType = "ENSEMBL", toType = "ENTREZID", 
            OrgDb = org.Bt.eg.db)$ENTREZID


pmd <- data[data$ENSEMBL %in% genes,]
pmd[,3:ncol(pmd)] <- log10(pmd[,3:ncol(pmd)] + 1)
pmd[,3:ncol(pmd)][pmd[,3:ncol(pmd)] < 0] <- 0
pmd$adult <- apply(pmd[,5:8], 1, mean)
pmd$fetal <- apply(pmd[,3:4], 1, mean)

non_pmd <- data[data$ENSEMBL %in% genes3,]
non_pmd[,3:ncol(non_pmd)] <- log10(non_pmd[,3:ncol(non_pmd)] + 1)
non_pmd[,3:ncol(non_pmd)][non_pmd[,3:ncol(non_pmd)] < 0] <- 0
non_pmd$adult <- apply(non_pmd[,5:8], 1, mean)
non_pmd$fetal <- apply(non_pmd[,3:4], 1, mean)

pmd <- pmd[,2:ncol(pmd)]
inputs <- reshape2::melt(pmd, measure.vars = c("adult", "fetal"),
                         variable.name = "Sample",value.name = "x")
inputs$group <- "PMD"
ggplot(inputs, aes(Sample, x)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 5))

non_pmd <- non_pmd[,2:ncol(non_pmd)]
inputs2 <- reshape2::melt(non_pmd,
                measure.vars = c("adult", "fetal"),
                variable.name = "Sample",value.name = "x")
inputs2$group <- "non_PMD"

all <- rbind(inputs, inputs2) %>% 
  select(c("ENSEMBL", "Sample", "x", "group"))

pdf("./methylation/pmd/pmd_flank.rna.pdf", 4, 4)
ggplot(all, aes(factor(Sample, levels = c("fetal", "adult")), 
                x, fill=group)) +
  geom_boxplot(color = "black", size = 0.7, width = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x= "", y = expression(log[10](TPM+1))) +
  cowplot::theme_cowplot() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top",
        legend.title = element_text(colour="black", size=10, face="bold")) +
  scale_fill_manual(values = c("#344ea2", "#ed2224"),
                    labels = c("non-PMD", "PMD"))
dev.off()

other <- data[!data$ENSEMBL %in% c(genes, genes3),]
other[,3:ncol(other)] <- limma::normalizeBetweenArrays(other[,3:ncol(other)],
)
other$fetal <- log10(apply(other[,3:4], 1, mean) + 1)
other$adult <- log10(apply(other[,5:8], 1, mean) + 1)

pdf("genes_not_involved.pdf", 4, 4)
other %>% 
  select(c("fetal", "adult")) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(width = 0.5, size = 0.6) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x= "", y = expression(log[10](TPM+1))) +
  cowplot::theme_cowplot() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top",
        legend.title = element_text(colour="black", size=10, face="bold")) +
  scale_fill_manual(values = c("#344ea2", "#ed2224"),
                    labels = c("non-PMD", "PMD"))
dev.off()

