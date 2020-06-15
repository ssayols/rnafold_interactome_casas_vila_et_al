#####################################
##
## Answer reviewer questions regarding specific biophysical properties of the 
## folds, that could attract specific proteins to bind:
##
##  * Rank top proteins binding to 5’, 3’ and CDR loops (by number of folds they bind to).
##    Are they overlapping or specific to the fold location?
##  * Are the fold location (5’, 3’ and CDR) having different biophysical properties?
##    G+C content, MFE from RNAfold for each structure classified by 5’, 3’ and UTR.
##    Is there any biophysical feature that contributes to this separation?
##  * How specific are these in vitro interactions? Or are these simply some high abundant
##    yeast proteins that somehow have a high affinity for folded RNA? A plot with
##    RBP abundance vs. number of folds it binds to
##
#####################################
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(ShortRead)
library(UpSetR)

CORES=16
set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/construct_biases.pdf", onefile=TRUE)

##
##  * Rank top proteins binding to 5’, 3’ and CDR loops (by number of folds they bind to).
##    Are they overlapping or specific to the fold location?
##

## read & dedup binders
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv")    # our list of binders
binders <- by(binders, 1:nrow(binders), function(x) {
    data.frame(X=unlist(strsplit(x$X, ";")),
               x[,!grepl("^(X|name)", colnames(x))],    # gene name is no longer valid for XXX;YYY rows
               row.names=NULL)
})
binders <- do.call(rbind, binders)
rownames(binders) <- binders[, 1]
binders <- as.matrix(binders[, -1])

# annotate binders
ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters=c("ensembl_gene_id"), values=rownames(binders),
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))

# read constructs names
constructs <- read.csv("constructs/constructs.csv")

# count, for each RBP, how many 3', 5' and CDS it binds to
p <- constructs$pos[match(colnames(binders), constructs$name)]
binders_counts_by_pos <- apply(binders, 1, function(x) factor(p[x > 0], levels=names(table(constructs$pos))))
binders_counts_by_pos <- do.call(rbind, lapply(binders_counts_by_pos, table))

# Venns with all and top10 binders of each type
upset(fromList(list(CDS =rownames(binders_counts_by_pos)[rev(order(binders_counts_by_pos[, "CDS" ]))[1:10]],
                    UTR3=rownames(binders_counts_by_pos)[rev(order(binders_counts_by_pos[, "UTR3"]))[1:10]],
                    UTR5=rownames(binders_counts_by_pos)[rev(order(binders_counts_by_pos[, "UTR5"]))[1:10]])))

# calculate the counts normalized by the number of constructs, and annotate with common gene name
binders_counts_by_pos <- cbind(as.data.frame(binders_counts_by_pos),
                               t(apply(binders_counts_by_pos, 1, function(x) x / table(constructs$pos))),
                               gene_name=ann$external_gene_name[match(rownames(binders_counts_by_pos), ann$ensembl_gene_id)])
colnames(binders_counts_by_pos)[4:6] <- paste0(colnames(binders_counts_by_pos)[4:6], "_norm")

##
##  * Are the fold location (5’, 3’ and CDR) having different biophysical properties?
##    G+C content, MFE from RNAfold for each structure classified by 5’, 3’ and UTR.
##    Is there any biophysical feature that contributes to this separation?
##

# Read fold sequences and calculate GC and MFE of the fold
x <- ShortRead::readFasta("./primers/feat.fasta")
constructs$seqs <- as.character(sread(x))[match(paste(constructs$name1, constructs$strand, constructs$name2),
                                                gsub("^ ", "", id(x)))]
constructs$GC  <- sapply(constructs$seq, function(x) seqinr::GC(unlist(strsplit(x, ""))))
constructs$MFE <- unlist(mclapply(constructs$seq, function(x) {
    rnafold <- system(sprintf("echo %s | bin/ViennaRNA/bin/RNAfold", x), intern=TRUE)
    as.numeric(gsub(".+\\((\\-*\\d+\\.\\d+)\\)$", "\\1", rnafold[2]))
}, mc.cores=CORES))

# do boxplots to see if there's differences between fold locations and GC/MFE
x <- reshape::melt(constructs, id.vars="pos", measure.vars=c("GC", "MFE"))
ggplot(x, aes(x=pos, y=value)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(comparisons=list(c("CDS", "UTR3"), c("CDS", "UTR5"), c("UTR5", "UTR3")), method="anova") +
    theme_bw() +
    facet_wrap(~variable, scales="free_y")

##
##  * How specific are these in vitro interactions? Or are these simply some high abundant
##    yeast proteins that somehow have a high affinity for folded RNA? A plot with
##    RBP abundance vs. number of folds it binds to
## 

# read and dedup bg
#bg <- filterIdentifiedProteins(filterWholeDataset(read.delim("./analysis/MQ/proteinGroups_reference_lysate.txt")))
bg <- read.delim("./analysis/MQ/proteinGroups_reference_lysate.rMQanalysis.txt")
bg <- do.call(rbind, by(bg, 1:nrow(bg), function(x) {
    data.frame(X=unlist(strsplit(x$Majority.protein.IDs, ";")),
               x[,colnames(x) != "Majority.protein.IDs"],
               row.names=NULL)
}))

# add abundance info to binders, and plot it vs. number of folds bound
binders_counts_by_pos$gene_id <- rownames(binders_counts_by_pos)
x <- reshape::melt(binders_counts_by_pos, id.vars="gene_id", measure.vars=c("CDS", "UTR3", "UTR5"))
x$abundance <- log10(bg$Intensity[match(x$gene_id, bg$X)])
ggplot(x, aes(x=value, y=abundance)) +
    geom_point() +
    geom_smooth(method="lm") +
    xlab("Binds to RNA folds") + ylab("Protein abundance (log10 intensity)")+
    theme_bw() +
    facet_wrap(~variable, scales="free_x")

write.csv(binders_counts_by_pos, file="analysis/construct_biases.csv")
dev.off()
