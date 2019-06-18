##
## plot the result of the GO term slim analysis from SGD
##
options(stringsAsFactors = FALSE)
library(ggplot2)
library(RColorBrewer)
library(reshape)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/GO")
pdf("GOterms.pdf")

# read SGD GO slim analysis results
go <- lapply(list(binders="binders vs all genes SGDslimmapper.txt", 
                  binders_no_rib="binders noRIB vs all genes SGDslimmapper.txt", 
                  all_rbp="RBPs vs all genes SGDslimmapper.txt", 
                  lysate_minus_rbp="lysate minus RBPs vs all genes SGDslimmaper.txt"), 
             function(f) {
    x <- read.delim(f)

    x$GO.term       <- gsub("(^\\s+|\\s+$)", "", x$GO.term)
    x$observed      <- as.numeric(gsub("(\\d+) out of (\\d+) genes, .+", "\\1", x$Frequency))
    x$totalObserved <- as.numeric(gsub("(\\d+) out of (\\d+) genes, .+", "\\2", x$Frequency))
    x$expected      <- as.numeric(gsub("(\\d+) of (\\d+) genes, .+", "\\1", x$Genome.Frequency))
    x$totalExpected <- as.numeric(gsub("(\\d+) of (\\d+) genes, .+", "\\2", x$Genome.Frequency))
    x$ratio         <- log2(x$observed / x$totalObserved) - log2(x$expected / x$totalExpected)

    filt <- grepl("other", x$GO.term) | x$expected == 0
    x[!filt, c("GO.term", "ratio")]
})

x <- Reduce(function(x, y) merge(x, y, by=1, all=T), go)
colnames(x) <- c("GO.term", names(go))

# plot as a barplot + heatmap
doBarplot <- function(x, i=10) {
    top <- apply(x[, -1], 2, function(x) rev(order(x)))
    i <- unique(as.integer(top[1:i, ]))
    df <- melt.data.frame(x[i, ])
    df$variable <- factor(df$variable, levels=rev(colnames(x)[-1]))
    df$value    <- ifelse(is.infinite(df$value) & df$value < 0, floor(min(df$value[!is.infinite(df$value)])), df$value)
    df$value    <- ifelse(is.infinite(df$value) & df$value > 0, ceiling(max(df$value[!is.infinite(df$value)])), df$value)
    df$GO.term  <- factor(df$GO.term, levels=df$GO.term[order(df$value[df$variable == colnames(x)[2]])])

    # barplot
    p <- ggplot(df, aes(x=GO.term, y=value, fill=variable)) + 
         geom_bar(stat="identity", position="dodge") +
         scale_fill_brewer(palette="Set2") +
         scale_y_continuous("log2 O/E ratio") +
         scale_x_discrete("") +
         coord_flip() + theme_bw() +
         theme(text=element_text(size=7), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7))
    print(p)

    # heatmap
    m <- as.matrix(x[i, -1])
    m <- ifelse(is.infinite(m) & m < 0, floor(min(m[!is.infinite(m)])), m)
    m <- ifelse(is.infinite(m) & m > 0, ceiling(max(m[!is.infinite(m)])), m)
#    o <- order(x[i, 1])
    o <- order(m[, 1])
    heatmap.2(m[o,], trace="none", scale="none", labRow=x[i, 1][o], cexCol=.6, cexRow=.6, notecex=.6, Rowv=NA, Colv=NA,
              dendrogram="none", col=colorRampPalette(brewer.pal(9, "Blues"))(100), cellnote=round(m[o,], 2),
              notecol="black", density.info="none", colsep=1:ncol(m), rowsep=1:nrow(m), sepwidth=c(0.05,0.05),
              key=FALSE, margin=c(5, 15), lmat=rbind(4:3, 2:1), lwid=c(1, 99), lhei=c(1,99))
}

doBarplot(x[, c("GO.term", "binders", "all_rbp", "lysate_minus_rbp")])  # all binders
doBarplot(x[, c("GO.term", "binders_no_rib", "all_rbp", "lysate_minus_rbp")])  # all binders but ribosomal

dev.off()
