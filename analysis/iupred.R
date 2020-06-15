#################################
##
## Calculate if our RBP are biases towards disordered regions.
##    1-get all proteins' sequence from biomart
##    2-call iupred and calculate % of bases with prob of disorder > 50% (as per pfam)
##    3-plot 3 lines with % of disordered bases: predicted (all proteome), measured proteome, identified (RBP)
##
#################################
library(biomaRt)
library(parallel)
library(rMQanalysis)
library(ggplot2)
library(reshape)

CORES=16
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/iupred.pdf", onefile=TRUE)

##
## 0-read input data
##
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv")    # our list of binders
bg      <- filterIdentifiedProteins(filterWholeDataset(read.delim("./analysis/MQ/proteinGroups_reference_lysate.txt")))

# supress proteins binding to COX17 control from binders list
binders <- binders[!binders$name %in% c('GCD11', 'SUI3', 'LSG1', 'PUF3'), ]

# deduplicate binders
binders <- by(binders, 1:nrow(binders), function(x) {
    data.frame(X=unlist(strsplit(x$X, ";")),
               x[,!grepl("^(X|name)", colnames(x))],    # gene name is no longer valid for XXX;YYY rows
               row.names=NULL)
})
binders <- do.call(rbind, binders)

# deduplicate bg
bg <- by(bg, 1:nrow(bg), function(x) {
    data.frame(X=unlist(strsplit(x$Majority.protein.IDs, ";")),
               x[,colnames(x) != "Majority.protein.IDs"],
               row.names=NULL)
})
bg <- do.call(rbind, bg)
bg <- bg[!(bg$X %in% binders$X), ]  # remove RBPs from the bg

##
##    1-get all proteins' sequence from biomart
##
seq <- getBM(attributes=c("ensembl_gene_id", "peptide"),
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
seq <- unique(seq, MARGIN=1)
seq <- seq[seq$peptide != "Sequence unavailable", ]
seq <- seq[sapply(seq$peptide, nchar) > 70, ]  # remove proteins shorter than 70aa (~250)

##
##    2-call iupred and calculate % of bases with prob of disorder > 50% (as per pfam)
##
iupred <- do.call(rbind, mclapply(setNames(seq$peptide, nm=seq$ensembl_gene_id), function(x) {
  writeLines(x, f <- tempfile())
  x <- read.table(text=system(paste("python3 ./bin/iupred2a/iupred2a.py", f, "long"), intern=TRUE), header=FALSE, sep="\t")
  sum(x$V3 > .5) / nrow(x)  # percentage of disordered bases
}, mc.cores=CORES))

##
##    3-plot 3 lines with % of disordered bases: predicted (all proteome), measured proteome, identified (RBP)
##
iupred <- data.frame(protein =rownames(iupred),
                     disorder=iupred,
                     class1  =ifelse(rownames(iupred) %in% binders$X, "identified",
                              ifelse(rownames(iupred) %in% bg$X     , "measured", "not_measured")),
                     class2  =ifelse(rownames(iupred) %in% binders$X, "identified", "predicted")
                     )
par(mfrow=c(2, 2))
boxplot(iupred$disorder ~ iupred$class1, cex=.6,
        main=paste("pvalue=", format.pval(anova(lm(iupred$disorder ~ iupred$class1))$`Pr(>F)`[1])))
boxplot(iupred$disorder ~ iupred$class2, cex=.6,
        main=paste("pvalue=", format.pval(anova(lm(iupred$disorder ~ iupred$class2))$`Pr(>F)`[1])))

# summarize (average) percentage of disorder in 100 bins
df1 <- lapply(split(iupred, iupred$class1), function(x) {
  x <- x[order(x$disorder), ]
  x$bin <- as.numeric(cut(1:nrow(x), 100))
  tapply(x$disorder, x$bin, mean)
})
df1 <- melt(df1)
df2 <- lapply(split(iupred, iupred$class2), function(x) {
  x <- x[order(x$disorder), ]
  x$bin <- as.numeric(cut(1:nrow(x), 100))
  tapply(x$disorder, x$bin, mean)
})
df2 <- melt(df2)
#ggplot(df, aes(x=indices, y=value, color=L1)) +
#  geom_point() +
ggplot(mapping=aes(x=indices, y=value, color=L1)) +
  geom_point(data=subset(df1, L1 %in% c("identified", "measured"))) +
  geom_point(data=subset(df2, L1 %in% c("predicted"))) +
  scale_color_discrete("") +
  xlab("") + ylab("Percentage of disorder") +
  theme_bw()

dev.off()
