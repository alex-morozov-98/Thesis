if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(enrichplot)
library(ggplot2)

### Select differentially expressed genes and draw volcano plot
# !NB: if you have previously analysed the file using Microsoft Excel and saved it, please, use sep = ';'
# Otherwise you should use sep = ','
foldc <- read.csv(file = 'C:/Users/Александр/Desktop/Data/R/final/local_fold_0vs3.csv', sep = ';')
head(foldc)
foldc$diffexpressed[foldc$avg_log2FC > 1.5 & foldc$p_val < 0.05] <- "UP"
foldc$diffexpressed[foldc$avg_log2FC < -1.5 & foldc$p_val < 0.05] <- "DOWN"
foldc$diffexpressed[is.na(foldc$diffexpressed)] <- "NO"
head(foldc)

p <- ggplot(data=foldc, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed)) + 
  geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-1.5, 1.5), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")
mycolors <- c("blue", "red", "gray")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
ggsave(filename='volcano_plot_0vs3.png',
       path = 'C:/Users/Александр/Desktop/Data/R/final/',
       plot = p3,
       width = 1000, 
       height = 690,
       limitsize = FALSE,
       scale = 3,
       units = "px")

up = subset(foldc, subset = diffexpressed == "UP")
down = subset(foldc, subset = diffexpressed == "DOWN")
write.table(up)
write.csv(up,'C:/Users/Александр/Desktop/Data/R/final/upregulated_0vs3.csv', row.names = TRUE)
write.table(down)
write.csv(down,'C:/Users/Александр/Desktop/Data/R/final/downregulated_0vs3.csv', row.names = TRUE)


### Preprocession of differentially expressed genes
# Upregulated
#setup <- read.csv(file = 'C:/Users/Александр/Desktop/Data/R/final/upregulated_0vs3.csv', sep=";" )
setup <- up
setup$FeatureName
entrezIDsup <- bitr(setup$FeatureName, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
head(entrezIDsup)
length(unique(entrezIDsup$SYMBOL))
length(unique(entrezIDsup$ENTREZID))
# Downregulated
#setdown <- read.csv(file = 'C:/Users/Александр/Desktop/Data/R/final/downregulated_0vs3.csv', sep=";" )
setdown <- down
setdown$FeatureName
entrezIDsdown <- bitr(setdown$FeatureName, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
head(entrezIDsdown)
length(unique(entrezIDsdown$SYMBOL))
length(unique(entrezIDsdown$ENTREZID))

### GO enrichment
# Upregulated 
EGO = enrichGO(gene = entrezIDsup$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 10,
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
EGO.df=as.data.frame(EGO@result[EGO@result$p.adjust <= 0.05,])
EGO_GO_up <- EGO
write.table(EGO_GO_up)
write.csv(EGO_GO_up,'C:/Users/Александр/Desktop/Data/R/final/EGO_GO_up.csv', row.names = TRUE)
# Downregulated
EGO = enrichGO(gene = entrezIDsdown$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 10,
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
EGO.df=as.data.frame(EGO@result[EGO@result$p.adjust <= 0.05,])
EGO_GO_down <- EGO
write.table(EGO_GO_down)
write.csv(EGO_GO_down,'C:/Users/Александр/Desktop/Data/R/final/EGO_GO_down.csv', row.names = TRUE)

### KEGG enrichment
# Upregulated
EGO_KEGG=enrichKEGG(entrezIDsup$ENTREZID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", 
                   keyType = 'kegg',qvalueCutoff=0.05)
EGO_KEGG=setReadable(EGO_KEGG, 'org.Hs.eg.db', 'ENTREZID')
EGO_KEGG.df=as.data.frame(EGO_KEGG@result[EGO_KEGG@result$p.adjust <= 0.05,])
EGO_KEGG_up <- EGO_KEGG
write.table(EGO_KEGG_up)
write.csv(EGO_KEGG_up,'C:/Users/Александр/Desktop/Data/R/final/EGO_KEGG_up.csv', row.names = TRUE)
# Downregulated 
EGO_KEGG=enrichKEGG(entrezIDsdown$ENTREZID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", 
                    keyType = 'kegg',qvalueCutoff=0.05)
EGO_KEGG=setReadable(EGO_KEGG, 'org.Hs.eg.db', 'ENTREZID')
EGO_KEGG.df=as.data.frame(EGO_KEGG@result[EGO_KEGG@result$p.adjust <= 0.05,])
EGO_KEGG_down <- EGO_KEGG
write.table(EGO_KEGG_down)
write.csv(EGO_KEGG_down,'C:/Users/Александр/Desktop/Data/R/final/EGO_KEGG_down.csv', row.names = TRUE)

### Reactome pathway enrichment
# Upregulated
Reactome <- enrichPathway(gene=entrezIDsup$ENTREZID,pvalueCutoff=0.05, readable=T)
Reactome.df=as.data.frame(Reactome@result[Reactome@result$p.adjust <= 0.05,])
Reactome_up <- Reactome
write.table(Reactome_up)
write.csv(Reactome_up,'C:/Users/Александр/Desktop/Data/R/final/Reactome_up.csv', row.names = TRUE)
# Downregulated
Reactome <- enrichPathway(gene=entrezIDsdown$ENTREZID,pvalueCutoff=0.05, readable=T)
Reactome.df=as.data.frame(Reactome@result[Reactome@result$p.adjust <= 0.05,])
Reactome_down <- Reactome
write.table(Reactome_down)
write.csv(Reactome_down,'C:/Users/Александр/Desktop/Data/R/final/Reactome_down.csv', row.names = TRUE)

### DOSE enrichment
# Upregulated
EGO_DO=enrichDO(entrezIDsup$ENTREZID, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,readable = T)
EGO_DO.df=as.data.frame(EGO_DO@result[EGO_DO@result$p.adjust <= 0.05,])
EGO_DO_up <- EGO_DO
write.table(EGO_DO_up)
write.csv(EGO_DO_up,'C:/Users/Александр/Desktop/Data/R/final/EGO_DO_up.csv', row.names = TRUE)
# Downregulated
EGO_DO=enrichDO(entrezIDsdown$ENTREZID, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,readable = T)
EGO_DO.df=as.data.frame(EGO_DO@result[EGO_DO@result$p.adjust <= 0.05,])
EGO_DO_down <- EGO_DO
write.table(EGO_DO_down)
write.csv(EGO_DO_down,'C:/Users/Александр/Desktop/Data/R/final/EGO_DO_down.csv', row.names = TRUE)

### visualize the results of analysis
types <- c("EGO_GO_up",'EGO_GO_down', 'EGO_KEGG_up', 'EGO_KEGG_down', 'EGO_DO_up', 'EGO_DO_down', 'Reactome_up', 'Reactome_down')
for (i in types)
{
  tryCatch({
    
    enrich <- eval(parse(text = i))

    cnet <- cnetplot(enrich, node_label="all", showCategory=20)
    file_name <- paste(i, "cnetplot", sep="_")
    ggsave(filename=paste(file_name, "png", sep="."),
       path = 'C:/Users/Александр/Desktop/Data/R/final/enrichment_pictures',
       plot = cnet,
       width = 1000, 
       height = 690,
       limitsize = FALSE,
       scale = 3,
       units = "px")

    dot <- dotplot(enrich, x="count", showCategory=20)
    file_name <- paste(i, "dotplot", sep="_")
    ggsave(filename=paste(file_name, "png", sep="."),
       path = 'C:/Users/Александр/Desktop/Data/R/final/enrichment_pictures',
       plot = dot,
       width = 1000, 
       height = 690,
       limitsize = FALSE,
       scale = 3,
       units = "px")

    heat <- heatplot(enrich, showCategory=20)
    file_name <- paste(i, "heatplot", sep="_")
    ggsave(filename=paste(file_name, "png", sep="."),
       path = 'C:/Users/Александр/Desktop/Data/R/final/enrichment_pictures',
       plot = heat,
       width = 1000, 
       height = 690,
       limitsize = FALSE,
       scale = 3,
       units = "px")
    
    edox2 <- pairwise_termsim(enrich)
    tree <- treeplot(edox2)
    file_name <- paste(i, "treeplot", sep="_")
    ggsave(filename=paste(file_name, "png", sep="."),
           path = 'C:/Users/Александр/Desktop/Data/R/final/enrichment_pictures',
           plot = tree,
           width = 1000, 
           height = 690,
           limitsize = FALSE,
           scale = 3,
           units = "px")

    edo <- pairwise_termsim(enrich)
    enrich_map <- emapplot(edo, showCategory=20)
    file_name <- paste(i, "enrichment_map", sep="_")
    ggsave(filename=paste(file_name, "png", sep="."),
       path = 'C:/Users/Александр/Desktop/Data/R/final/enrichment_pictures',
       plot = enrich_map,
       width = 1000, 
       height = 690,
       limitsize = FALSE,
       scale = 3,
       units = "px")}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


