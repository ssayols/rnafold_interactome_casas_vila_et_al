#############################################
##
## Take the results from heatmap and redo the plot,
## this time naming the rows and the columns witht the
## protein name instead of id
##
#############################################
options(stringsAsFactors=F)
library(gplots)
library(RColorBrewer)
library(biomaRt)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/heatmap2.pdf", onefile=TRUE, width=14, height=14)

# read cosntructs info
constructs <- read.csv("./constructs/constructs.csv")
constructs$name <- paste0(constructs$plate, constructs$position)
constructs$pos  <- gsub(".+_(.+)_.$", "\\1", constructs$name1)

ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters=c("ensembl_gene_id"), values=constructs$gene,
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
constructs$gene_name <- ann$external_gene_name[match(constructs$gene, ann$ensembl_gene_id)]
constructs$gene_name <- ifelse(constructs$gene_name == "", constructs$gene, constructs$gene_name)

# read heatmap.R results
m <- read.csv("./analysis/heatmap_1+1_threshold.csv")
m_row_names <- m$name
m <- as.matrix(m[, grepl("^P", colnames(m))])
m_col_names <- constructs$gene_name[match(colnames(m), constructs$name)]

unique.binders <- apply(m, 1, function(x) sum(x > 0) == 1)

# do plot
pal <- c(colorRampPalette(c(rep("#E41A1C", 2), "black"))(50),    # controls, in sort of read
         colorRampPalette(c("black", rep("yellow", 2)))(50))     # binders, in yellow
pos <- factor(constructs$pos[match(colnames(m), constructs$name)], levels=c("UTR5", "CDS", "UTR3"))
cc  <- c("#4DAF4A", "#E41A1C", "#377EB8")[pos]  # taken from the brewer.pal(3, "Set1") palette, but in different order
o <- order(pos, colSums(m))

h <- heatmap.2(m[, o], ColSideColors=cc[o], labCol=m_col_names[o], labRow=m_row_names, 
               RowSideColors=brewer.pal(3, "Set1")[as.factor(unique.binders)],
               Rowv=T, Colv=F, col=pal, dendrogram="none", cexRow=0.5, cexCol=0.4, key=T, srtCol=90,
               trace="none", density.info="none",
               main="Proteins with higher SILAC ratios than 1")

# save the heatmap matrix as it is shown
out <- as.data.frame(m[rev(h$rowInd), o])
out <- cbind(out, m_row_names[rev(h$rowInd)])  # heatmap reverses rows after clustering
out <- rbind(out, c(m_col_names[o], ""))       # we didn't cluster cols
write.csv(out, file="./analysis/heatmap2.csv", row.names=FALSE)

dev.off()
