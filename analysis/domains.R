##########################
##
## Analysis of binding domains from Tuschel's Nat Rev Gen
##
##   1-read input data
##   2-annotate proteins with the domain names from biomart
##   3-Tuschel's Pfam vs. GO tabulated plot of the binders
##   4-barplot summary
##   5-rescue the non RBD with human orthologs
##   6-overrepresentation test of Tuschel's pfam domains present in the binders
##   6B-similar representation, but counting in a barplot the ocurrences of Tuschel's domains
##   7-identify intensity of binders in the whole proteome reference lysate
##
##########################
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(rMQanalysis)

set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/domains.pdf", onefile=TRUE)

##
## 1-read input data
##
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv")    # our list of binders
bg      <- filterIdentifiedProteins(filterWholeDataset(read.delim("./analysis/MQ/proteinGroups_reference_lysate.txt")))
pfam    <- read.delim("./analysis/db/Pfam-A.clans.tsv", head=F)    # translation between ensembl pfam_id and Tuschel's pfam name
rbd.pfam  <- read.delim("./analysis/db/PfamA_tuschel.csv")
rbd.go    <- read.delim("./analysis/db/GO_tuschel.csv")

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

##
## 2-annotate proteins with the domain names from biomart
##
annotate.domains <- function(x, dataset) {

    # get gene's pfam domains and GO terms from biomart
    ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "pfam", "name_1006"),
                      filters=c("ensembl_gene_id"), values=x$X,
                      mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset=dataset))

    # check which x are RNA related from Tuschel's pfamA list
    x$pfam <- sapply(x$X, function(x) {
        pfam_id   <- ann$pfam[ann$ensembl_gene_id == x]
        pfam_name <- pfam$V4[match(pfam_id, pfam$V1, nomatch=0)]
        any(pfam_name %in% rbd.pfam$Pfam.RNA.binding.domains)
    })

    # check which x have the "RNA binding" GO term
    x$go_rna_bind <- sapply(x$X, function(x) any(grepl("\\bRNA binding\\b", ann$name_1006[ann$ensembl_gene_id == x])))

    # check which x are RNA related from Tuschel's GO list
    x$go_rna_rel <- sapply(x$X, function(x) any(ann$name_1006[ann$ensembl_gene_id == x] %in% rbd.go$selected.GO.terms.for.RN..related.processes))

    x
}

binders <- annotate.domains(binders, "scerevisiae_gene_ensembl")
bg      <- annotate.domains(bg     , "scerevisiae_gene_ensembl")

##
## 3-Tuschel's Pfam vs. GO tabulated plot of the binders
##
x <- data.frame(pfam=ifelse(binders$pfam, "known", "unknown"),
                go=ifelse(binders$go_rna_bind, "rna bind", ifelse(binders$go_rna_rel, "rna related", "unrelated")))
x <- table(x$pfam, x$go)
#mosaicplot(x, xlab="Tuschel's PFAM census", ylab="Tuschel's GO census", main="", col=brewer.pal(3,"Set1"))
x <- melt(x)

ggplot(x, aes(x=Var2, y=Var1)) +
    geom_point(aes(size=value), color="red", alpha=.5) +
    geom_text(aes(label=value)) +
    scale_size(range = c(5, 50), guide=FALSE) +
    scale_x_discrete("Tuschel's GO census") +
    scale_y_discrete("Tuschel's PFAM census") +
    theme_minimal()

##
## 4-barplot summary
##
pal <- brewer.pal(5, "Set3")
summ <- function(x, what) {
    data.frame(both=sum(x$pfam & x$go_rna_bind),
               pfam=sum(x$pfam & !x$go_rna_bind),
               go_rna_bind=sum(!x$pfam & x$go_rna_bind),
               go_rna_related=sum(!x$pfam & !x$go_rna_bind & x$go_rna_rel),
               unrelated=sum(!x$pfam & !x$go_rna_bind & !x$go_rna_rel),
               what=what)
}
x <- rbind(summ(binders, "binders"),
           summ(bg, "bg"))
x <- melt(x)
x$variable <- factor(x$variable,
                     levels=c("pfam", "go_rna_bind", "both", "go_rna_related", "unrelated"),
                     labels=c("pfam", "GO:rna bind", "both", "GO:rna related", "unrelated"))
ggplot(x, aes(x=what, y=value, fill=variable, label=value)) +
    geom_bar(stat="identity", position="fill") + #, alpha=.5) +
    geom_text(position="fill", vjust="inward") +
    theme_minimal() +
    scale_fill_manual(values=pal[1:5]) +       # unrelated always pal[5]
    scale_x_discrete("") + scale_y_continuous("")

# and now, divide into PFAM and GO annotated (2 plots)
x$pfam <- factor(ifelse(as.character(x$variable) %in% c("pfam", "both"), "pfam", "unrelated"),
                 levels=c("pfam", "unrelated"),
                 labels=c("RBD", "no RBD"))
x2 <- do.call(rbind, by(x, x[, c("what", "pfam")], function(x) {
  res <- x[1, , drop=F]
  res$value <- sum(x$value)
  res
}))
ggplot(x2, aes(x=what, y=value, fill=pfam, label=value)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(position="fill", vjust="inward") +
    theme_minimal() +
    scale_fill_manual(values=pal[c(1, 5)]) +   # unrelated always pal[5]
    scale_x_discrete("") + scale_y_continuous("")

x$GO <- factor(ifelse(as.character(x$variable) %in% c("GO:rna bind", "both"), "bind", 
               ifelse(as.character(x$variable) %in% c("GO:rna related"), "related", "unrelated")),
               levels=c("bind", "related", "unrelated"),
               labels=c("RNA binding", "RNA related", "unrelated"))
x2 <- do.call(rbind, by(x, x[, c("what", "GO")], function(x) {
  res <- x[1, , drop=F]
  res$value <- sum(x$value)
  res
}))
ggplot(x2, aes(x=what, y=value, fill=GO, label=value)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(position="fill", vjust="inward") +
    theme_minimal() +
    scale_fill_manual(values=pal[c(2, 4, 5)]) +   # unrelated always pal[5]
    scale_x_discrete("") + scale_y_continuous("")

##
## 5-rescue the non RBD with human orthologs
##
# get human orthologs of the binders with non RNA related domains
binders.unrelated <- binders$X[!binders$pfam & !binders$go_rna_bind & !binders$go_rna_rel]
orth <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene",
                           "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_orthology_confidence"),
              filters=c("ensembl_gene_id"), values=binders.unrelated,
              mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))

# annotate the human domains
colnames(orth)[grepl("hsapiens_homolog_ensembl_gene", colnames(orth))] <- "X"
orth <- annotate.domains(orth, "hsapiens_gene_ensembl")

# collapse rows by saccer ensembl gene id
orth <- by(orth, orth$ensembl_gene_id, function(x) {
    data.frame(ensembl_gene_id=paste(unique(x$ensembl_gene_id), collapse=";"),
               external_gene_name=paste(unique(x$external_gene_name), collapse=";"),
               hsapiens_homolog_ensembl_gene=paste(unique(x$X), collapse=";"),
               hsapiens_homolog_associated_gene_name=paste(unique(x$hsapiens_homolog_associated_gene_name), collapse=";"),
               pfam=any(x$pfam),
               go_rna_bind=any(x$go_rna_bind),
               go_rna_rel=any(x$go_rna_rel))
})
orth <- do.call(rbind, orth)

# add the binders.unrelated which have no human homolog (SINCE WORKING WITH GENE ID, SHOULD NOT NECESSARY ANYMORE)
orth <- merge(as.data.frame(binders.unrelated), orth, by.x=1, by.y=0, all.x=TRUE)
orth$pfam <- ifelse(is.na(orth$pfam), FALSE, orth$pfam)
orth$go_rna_bind <- ifelse(is.na(orth$go_rna_bind), FALSE, orth$go_rna_bind)
orth$go_rna_rel  <- ifelse(is.na(orth$go_rna_rel ), FALSE, orth$go_rna_rel )

# barplot summary
x <- melt(summ(orth, "human orthologs"))
x$variable <- factor(x$variable, levels=c("pfam", "go_rna_bind", "both", "go_rna_related", "unrelated"),
                     labels=c("pfam", "GO:rna bind", "both", "GO:rna related", "unrelated"))
x <- x[rev(order(x$variable)),]
ggplot(x, aes(x=what, y=value, fill=variable, label=value)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(position="fill", vjust="inward") +
    theme_minimal() +
    scale_fill_manual(values=pal[1: 5]) +   # unrelated always pal[5]
    scale_x_discrete("") + scale_y_continuous("")

# save results
write.csv(binders, file="./analysis/domains_binders.csv", row.names=F)
write.csv(orth, file="./analysis/domains_binders_unrelated.csv", row.names=F)

##
##   6-Overrepresentation test of Tuschel's pfam domains present in the binders
##
# get the whole genes as the bg
#bg2 <- getBM(attributes=c("ensembl_gene_id"), mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
#colnames(bg2) <- "X"

# get gene's pfam domains from biomart
#pfam.count <- lapply(list(binders$X, bg2$X), function(x) {
pfam.count <- lapply(list(binders$X, bg$X), function(x) {
    ann <- getBM(attributes=c("ensembl_gene_id", "pfam"), filters="ensembl_gene_id", values=x,
                 mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
    ann <- unique(ann, MARGIN=1)    # don't count multiple times proteins with 2+ an X domain

    # count number of domains from Tuschel's pfam census in the list of binders
    sapply(rbd.pfam$Pfam.RNA.binding.domains, function(x) {
        sum(ann$pfam %in% pfam$V1[pfam$V4 == x])
    })
})

# remove RBD not associated with any protein in yest
out <- pfam.count[[1]] == 0 & pfam.count[[2]] == 0
pfam.count <- lapply(pfam.count, function(x) x[!out])

# calculate the significance of the recovery
# WARNING: I realized the significance is heavily biased by the fact that sometimes proteins
#          where detected in the binders but not in the bg, which is a technical MS issue...
significance <- mapply(function(x, y, binders, bg) {
    phyper(x, y, bg, binders, lower.tail=FALSE)   # where x=binders+domain, y=bg+domain
}, pfam.count[[1]], pfam.count[[2]], MoreArgs=list(binders=sum(pfam.count[[1]]), bg=sum(pfam.count[[2]])), SIMPLIFY=F)
significance <- p.adjust(unlist(significance), method="BH")
significance <- ifelse(significance == 0, .001, significance)   # fix "infinite" enrichments

# plot something
x <- data.frame(domain=names(pfam.count[[1]]),
                binders=pfam.count[[1]],
                bg=pfam.count[[2]],
                significance=-log10(significance),
                recovery=100 * pfam.count[[1]] / pfam.count[[2]])
x$recovery <- ifelse(x$recovery > 100, 100, x$recovery)

ggplot(x, aes(x=bg, y=binders)) +
    geom_point(aes(color=significance, size=recovery), alpha=.5) +
    geom_text_repel(data=subset(x, (binders > 1 & significance > 1) | binders > 2), mapping=aes(label=domain)) +
    scale_size("% recovery", range=c(1, 10)) +
    scale_color_distiller("-log10 pval", palette="Spectral") +
    scale_x_continuous("Background") +
    scale_y_continuous("Binders") +
    ggtitle("RBD in Tuschel's PFAM census") +
    theme_minimal()

# same, but as a barplot. Only the significant domains (pval < .1)
y <- melt(x, measure.vars=c("binders", "bg"))
ggplot(subset(y, significance > 1), aes(x=domain, y=value, fill=variable)) +
    geom_bar(stat="identity", position="dodge") +
    ylab("counts") + xlab("domain") +
    ggtitle("Significantly enriched RBD in Tuschel's PFAM census (fdr < .1)") +
    coord_flip() +
    theme_minimal()

##
##   6B-similar representation, but counting in a barplot the ocurrences of Tuschel's domains
##
df <- data.frame(domain=names(pfam.count[[1]][pfam.count[[1]] > 0]),
                 counts=pfam.count[[1]][pfam.count[[1]] > 0])
df <- df[order(df$counts, df$domain), ]
df$domain <- factor(df$domain, levels=df$domain)

ggplot(df, aes(x=domain, y=counts, fill=factor(counts))) +
    geom_bar(stat="identity") +
    scale_fill_brewer("", palette="Blues") +
    ggtitle("RBD in Tuschel's PFAM census") +
    coord_flip() +
    theme_minimal()

ggplot(subset(df, counts > 1), aes(x=domain, y=counts, fill=factor(counts))) +
    geom_bar(stat="identity") +
    scale_fill_brewer("", palette="Blues") +
    ggtitle("RBD in Tuschel's PFAM census") +
    coord_flip() +
    theme_minimal()

##
##   7-identify intensity of binders in the whole proteome reference lysate
##
# calculate intensity of the protein in the background. It's the average of the H+L experiments
#bg_intensities <- log10(bg$Intensity.heavy + bg$Intensity.light + 1) # that used to be in an old protein groups
#bg_intensities <- bg_intensities[bg_intensities > 0]    # remove not detected proteins
bg_intensities <- log10(bg$Intensity)
pal <- brewer.pal(3, "Set1")[1:2] # red, blue
is.binder <- sapply(bg$X, function(x) unlist(strsplit(x, ";")) %in% binders$X)
bg_intensities_binders <- bg_intensities[is.binder]
bg_intensities_nobinders <-  sample(bg_intensities[!is.binder], sum(is.binder))
boxplot(list(binders=bg_intensities_binders, other=bg_intensities_nobinders), notch=T,
        col=pal, main="intensity of binders in the whole proteome", ylab="Intensity")

# density of intensities
plot(density(bg_intensities_binders, na.rm="T"), col=pal[1], lwd=2,
     main="distribution of intensities")
lines(density(sample(bg_intensities_nobinders), na.rm="T"), col=pal[2], lwd=2)    # take only few values
legend("topright", legend=c("binders", "other"), fill=pal)

dev.off()

# plot the proteins sorted by intensity, and identify where the binders are
cairo_pdf("./analysis/domains4x4.pdf", onefile=TRUE, width=4, height=4)
o <- order(bg_intensities)
plot(0, xlim=c(1, length(bg_intensities)), ylim=range(bg_intensities), ylab="Intensity",
     main="Intensity vs. protein abundance index")
points(1:length(bg_intensities), bg_intensities[o], pch=16,
       col=ifelse(is.binder[o], paste0(pal[1], "50"), pal[2]),
       cex=ifelse(is.binder[o], 2, 1))
legend("topleft", legend=c("binders", "other"), fill=pal)

dev.off()
