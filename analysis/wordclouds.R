################################
##
## clouds of GO slim terms of the proteins in the 2 dimensions of the heatmap (constructs and binders)
##
################################
library(wordcloud)
library(RColorBrewer)
library(WriteXLS)

WD <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib"
setwd(WD)

##
## read heatmap matrix with constructs and binders
##
consts.dict  <- read.csv("./constructs/constructs.csv")
ensembl.dict <- read.delim("./functional_follow_up/ensembl_dic.txt")
goslim       <- read.delim("./analysis/db/go_slim_mapping.tab", comment.char="#")
goslim       <- goslim[goslim$GOtype == "P", ]  # only BP
goslim       <- goslim[goslim$term != "biological_process", ]
h <- read.csv("./analysis/heatmap_1+1_threshold.csv")
h <- h[-1:-4, ] # remove controls
binders      <- h$X
h <- as.matrix(h[, -1:-2])
constructs   <- consts.dict$gene[match(colnames(h), consts.dict$name)]

##
## create wordclouds with the binders each construct
##
# get the terms
get.terms <- function(ratios, genes) {
    genes.enriched <- unlist(strsplit(genes[ratios > 0], ";"))
    genes.enriched.SGD <- ensembl.dict$SGD.Gene[match(genes.enriched, ensembl.dict$Ensembl.Gene.ID)]
    terms <- goslim[goslim$id %in% genes.enriched.SGD, c("gene", "term"), drop=F]
    if(nrow(terms) > 0) {
        terms <- by(terms, terms$term, function(x) data.frame(term=x$term[1], genes=paste(sort(unique(x$gene)), collapse="; ")))
        terms <- do.call(rbind, terms)
    }

    list(terms=terms$term,
         genes=terms$genes,
         numgenes=length(genes.enriched),
         numterms=sum(!duplicated(goslim$term[goslim$id %in% genes.enriched.SGD])))
}

# do the wordcloud
do.plot <- function(x, main) {
    if(x$numterms > 0) {
        main <- paste("\n", unique(ensembl.dict$Associated.Gene.Name[match(main, ensembl.dict$Ensembl.Gene.ID)]),
                      "(", x$numgenes, "genes,", x$numterms, "terms )")
        wordcloud(x$terms, random.order=F, min.freq=2, colors=brewer.pal(8, "Dark2"))
        title(main)
    }
}


pdf("./analysis/wordclouds_binders.pdf")
par(mfrow=c(2, 2))
mapply(do.plot, x <- apply(h, 1, get.terms, constructs), strsplit(binders, ";"))
dev.off()
out <- mapply(function(x, id, symbol) cbind(x[[1]], x[[2]], id, symbol),
              x, binders, ensembl.dict$Associated.Gene.Name[match(binders, ensembl.dict$Ensembl.Gene.ID)])
out <- out[sapply(out, ncol) == 4]
out <- as.data.frame(do.call(rbind, out))
WriteXLS("out", ExcelFileName="./analysis/wordclouds_binders.xls")

pdf("./analysis/wordclouds_constructs.pdf")
par(mfrow=c(2, 2))
mapply(do.plot, x <- apply(h, 2, get.terms, binders), constructs)
dev.off()
out <- mapply(function(x, id, symbol) cbind(x[[1]], x[[2]], id, symbol),
              x, constructs, ensembl.dict$Associated.Gene.Name[match(constructs, ensembl.dict$Ensembl.Gene.ID)])
out <- out[sapply(out, ncol) == 4]
out <- as.data.frame(do.call(rbind, out))
WriteXLS("out", ExcelFileName="./analysis/wordclouds_constructs.xls")

