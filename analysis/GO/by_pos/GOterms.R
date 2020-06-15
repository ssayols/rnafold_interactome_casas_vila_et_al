############################
##
## GO terms associated with binders, split by UTR5-CDS-UTR3
##
############################
options(stringsAsFactors=FALSE)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/GO/by_pos/")

# read constructs and binder info
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv")
constructs <- read.csv("./constructs/constructs.csv")
constructs <- constructs[constructs$name %in% colnames(binders), ]

ann <- data.frame(gene_id  =unlist(strsplit(binders$X   , ";")),
                  gene_name=unlist(strsplit(binders$name, ";")))


## call GO (will use SGD's Gene Ontology Slim Term Mapper)
#urls <- tapply(constructs$name, constructs$pos, function(x) {
#    binders.pos <- binders$name[apply(binders[, x], 1, function(x) any(x > 0))]
#    list(BP=paste0("https://www.yeastgenome.org/goSlimMapper?",
#                  "genes=", paste(binders.pos, collapse="+"), "&",
#                  "uploadFile=&",
#                  "slim_type=Yeast+GO-Slim%3A+process&",
#                  "submit=Submit+Form"),
#         MF=paste0("https://www.yeastgenome.org/goSlimMapper?",
#                  "genes=", paste(binders.pos, collapse="+"), "&",
#                  "uploadFile=&",
#                  "slim_type=Yeast+GO-Slim%3A+function&",
#                  "submit=Submit+Form"))
#})
#
#invisible({
#    Map(urls, names(urls), f=function(x, x.name) {
#        cat("\npaste this url in the browser and download the result as tab-delimited:\n\n", x$BP,
#            "\n\nand save as GOterms.", x.name, ".BP.txt\n", sep="")
#        cat("paste this url in the browser and download the result as tab-delimited:\n\n", x$MF,
#            "\n\nand save as GOterms.", x.name, ".MF.txt\n", sep="")
#    })
#})

# need to use gene_id instead of common name, since some are ambiguous and SGD's Slirm Mapper complains
rbp <- tapply(constructs$name, constructs$pos, function(x) {
    unlist(strsplit(binders$X[apply(binders[, x], 1, function(x) any(x > 0))], ";"))
})

# upload the file in https://www.yeastgenome.org/goSlimMapper and run BP and MF.
# save the results as tab-delimited as GOterms.<pos>.[BP|MF].txt
Map(rbp, names(rbp), f=function(x, x.name) {
    write.table(x, file=paste0(x.name, ".binders.txt"), row.names=FALSE, quote=FALSE)
})


pdf("GOterms.pdf")

sapply(c("BP", "MF"), function(what) {
    # read SGD GO slim analysis results
    go <- lapply(list(CDS =paste0("GOterms.CDS.", what, ".txt"),
                      UTR3=paste0("GOterms.UTR3.", what, ".txt"),
                      UTR5=paste0("GOterms.UTR5.", what, ".txt")),
                 function(f) {
        x <- read.delim(f)
        x$ratio <- log2(x$NUM_LIST_ANNOTATIONS / x$LIST_SIZE) - log2(x$TOTAL_NUM_ANNOTATIONS / x$POPULATION_SIZE)

        x$ANNOTATED_GENES2 <- sapply(x$ANNOTATED_GENES, function(genes) paste(ann$gene_name[match(unlist(strsplit(genes, ", ")), ann$gene_id)], collapse=", "))
        write.table(x, file=f, row.names=FALSE, quote=FALSE, sep="\t")

        filt <- grepl("other", x$TERM) | x$TOTAL_NUM_ANNOTATIONS < 5   # exclude small categories with less than 5 genes
        x[!filt, c("TERM", "ratio")]
    })

    x <- Reduce(function(x, y) merge(x, y, by=1, all=T), go)
    colnames(x) <- c("GO.term", names(go))

    # plot as a barplot + heatmap
    doBarplot <- function(x, i=25) {
        top <- apply(x[, -1], 2, function(x) rev(order(abs(x))))
        i <- unique(as.integer(top[1:i, ]))
        df <- melt.data.frame(x[i, ])
        df$variable <- factor(df$variable, levels=rev(colnames(x)[-1]))
        df$value    <- ifelse(is.infinite(df$value) & df$value < 0, floor(min(df$value[!is.infinite(df$value)])), df$value)
        df$value    <- ifelse(is.infinite(df$value) & df$value > 0, ceiling(max(df$value[!is.infinite(df$value)])), df$value)
        df$GO.term  <- factor(df$GO.term, levels=df$GO.term[order(df$value[df$variable == colnames(x)[2]])])

        # barplot
        p <- ggplot(df, aes(x=GO.term, y=value, fill=variable)) + 
             geom_bar(stat="identity", position="dodge") +
             scale_fill_brewer("binds to", palette="Set2") +
             scale_y_continuous("log2 O/E ratio") +
             scale_x_discrete("") +
             ggtitle(what) +
             coord_flip() + theme_bw() +
             theme(text=element_text(size=7), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7))
        print(p)
    }

    doBarplot(x)
})

dev.off()
