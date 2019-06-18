##########################
##
## Analysis of modifications in known/putative binding domains
##
##   1-read input data
##   2-mark the modifications whether they belong to binders, and mark the modifications
##     if they fall into an annotated RBD in pfam
##   3-load the background from the acetyl/formyl, do some plots of this contingency table:
##            +------+------+
##            |  RBD | !RBD |
##    +-------+------+------+
##    |  mark |      |      |
##    +-------+------+------+
##    | !mark |      |      |
##    +-------+------+------+
##   4-same for the binders
##   **NOTE: bg will be proteins quantified and with known domains (~10% out). The domains will be classified as RBD/DBD/none
##   **NOTE: bg will we ALL --> bg without substracting binders
##
##########################
library(ggplot2)
library(reshape)
library(biomaRt)
library(WriteXLS)
library(rMQanalysis)

set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./functional_follow_up/acetyl_formyl.pdf", onefile=TRUE)

##
## 1-read input data
##
pfam     <- read.delim("./analysis/db/Pfam-A.clans.tsv", head=F) # translation between ensembl pfam_id and Tuschel's pfam name
dbd.pfam <- read.delim("./analysis/db/PfamA_sergi_DBD.csv", comment.char="#") # Sergi's pfam dbd census
rbd.pfam <- read.delim("./analysis/db/PfamA_tuschel.csv")        # Tuschel's pfam rbd census
rbd.pfam$pfam <- pfam$V1[match(rbd.pfam$Pfam.RNA.binding.domains, pfam$V4)]
binders  <- read.csv("./analysis/heatmap_1+1_threshold.csv")     # our list of binders
binders  <- unique(unlist(strsplit(binders$X, ";")))
bg       <- filterIdentifiedProteins(filterWholeDataset(read.delim("./functional_follow_up/acetyl_formyl_data_bg.txt")))
bg       <- unlist(strsplit(bg$Majority.protein.IDs, ";"))

l <- list(acetyl="./functional_follow_up/acetyl_data_short.csv",
          formyl="./functional_follow_up/formyl_data_short.csv")
modif   <- mapply(function(f, type) {
    x <- read.csv(f)
    colnames(x)    <- c("Proteins", "Positions.within.proteins", "Number.of.modified.AA",
                        "Sequence.window", "Probabilities", "Score.diffs", "Position.in.peptide")
    x$modification <- type
    x
}, l, names(l), SIMPLIFY=FALSE)
modif <- do.call(rbind, modif)

##
##   2-mark the modifications whether they belong to binders, and mark the modifications
##     if they fall into an annotated RBD in pfam
##
# deduplicate rows
modif <- by(modif, modif$Proteins, function(x) {
             data.frame(Proteins=unlist(strsplit(x$Proteins, ";")),
                        Positions.within.proteins=as.integer(unlist(strsplit(x$Positions.within.proteins, ";"))),
                        x[, -1:-2])
         })
modif <- do.call(rbind, modif)

# get the protein domains from biomart
domains <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "uniprot_swissprot",
                              "pfam", "pfam_start", "pfam_end"),
#                 filter="ensembl_gene_id",
#                 values=unique(modif$Proteins),
                 mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
domains <- domains[!is.na(domains$pfam_start) & !is.na(domains$pfam_end), ] # ~2000 domains with unknown start/end. Discard
domains$is.rbd <- domains$pfam %in% rbd.pfam$pfam
domains$is.dbd <- domains$pfam %in% dbd.pfam$Accession

# check if binder has a modification
binders.modif <- sapply(binders, grep, modif$Proteins)
binders.modif <- binders.modif[sapply(binders.modif, function(x) length(x) > 0)]
modif$is.binder <- FALSE
modif$is.binder[unique(unlist(binders.modif))] <- TRUE

# check if bg has a modification
bg.modif <- sapply(bg, grep, modif$Proteins)
bg.modif <- bg.modif[sapply(bg.modif, function(x) length(x) > 0)]
modif$is.bg <- FALSE
modif$is.bg[unique(unlist(bg.modif))] <- TRUE

# annotate the modification with the domain name
modif <- by(modif, modif$Proteins, function(m) {
    d <- domains[domains$ensembl_gene_id %in% m$Proteins[1], ]
    m$domain <- sapply(m$Positions.within.proteins, function(x) {
        i <- d$pfam_start <= x & x <= d$pfam_end
        if(any(i))
            paste(unique(d$pfam[i]), collapse=";")  # should be only 1...
        else
            "no domain annotated here"
        })
    m
})
modif <- do.call(rbind, modif)

# annotate if the domain is RBD. Consider modif$domain contains only 1 domain (no paste(x, collapse=";")
modif$is.rbd <- modif$domain %in% rbd.pfam$pfam
modif$is.dbd <- modif$domain %in% dbd.pfam$Accession

# add some useful names
modif$uniprot   <- domains$uniprot_swissprot [match(modif$Proteins, domains$ensembl_gene_id)]
modif$gene_name <- domains$external_gene_name[match(modif$Proteins, domains$ensembl_gene_id)]

WriteXLS(modif, ExcelFileName="./functional_follow_up/acetyl_formyl.xls")

##
##   3-load the background from the acetyl/formyl, do some plots of this contingency table:
##            +------+------+
##            |  RBD | !RBD |
##    +-------+------+------+
##    |  mark |      |      |
##    +-------+------+------+
##    | !mark |      |      |
##    +-------+------+------+
##   4-same for the binders
##   **NOTE: bg will be proteins quantified and with known domains (~10% out). The domains will be classified as RBD/DBD/none
##   **NOTE: bg will we ALL --> bg without substracting binders
##
bg <- bg[bg %in% domains$ensembl_gene_id]   # all proteins detected in the masspec, with a known domain in pfam

# summarize & barplot
summ <- function(x, what) {
    bd <- domains$ensembl_gene_id[domains$is.rbd | domains$is.dbd]  # genes with a known nucleotide binding domain
    data.frame(bd_nomodif      =sum(x %in% bd & !(x %in% modif$Proteins)),
               bd_modif_dbd    =sum(x %in% bd &  (x %in% modif$Proteins) &
                                    sapply(x, function(x)  any(domains$pfam[domains$ensembl_gene_id == x] %in%
                                                               modif$domain[modif$Proteins == x & ( modif$is.rbd |  modif$is.dbd)]))),
               bd_modif_nodbd  =sum(x %in% bd &  (x %in% modif$Proteins) &
                                    sapply(x, function(x) !any(domains$pfam[domains$ensembl_gene_id == x] %in%
                                                               modif$domain[modif$Proteins == x & ( modif$is.rbd | modif$is.dbd)]))),
               nobd_nomodif    =sum(!(x %in% bd) & !(x %in% modif$Proteins)),
               nobd_modif_nodbd=sum(!(x %in% bd) &  (x %in% modif$Proteins) &
                                    (sapply(x, function(x)  any(domains$pfam[domains$ensembl_gene_id == x] %in%
                                                                modif$domain[modif$Proteins == x & (!modif$is.rbd & !modif$is.dbd)])) |
                                     sapply(x, function(x) !any(domains$pfam[domains$ensembl_gene_id == x] %in%
                                                                modif$domain[modif$Proteins == x & (!modif$is.rbd & !modif$is.dbd)])))),
               what=what)
}
x <- rbind(summ(binders, "binders"),
           summ(bg, "bg"))
x <- melt(x)
x$variable <- factor(x$variable,
                     levels=c("bd_nomodif", "bd_modif_dbd", "bd_modif_nodbd", "nobd_modif_nodbd", "nobd_nomodif"),
                     labels=c("(+)DBD/RBD (-)modif", "(+)DBD/RBD (+)modif", "(+)DBD/RBD (-)modif, (-)DBD/RBD (+)modif",
                              "(-)DBD/RBD (+)modif", "(-)DBD/RBD (-)modif"))
x <- x[rev(order(x$variable)),]
ggplot(x, aes(x=what, y=value, fill=variable, label=value)) +
    geom_bar(stat="identity", position="fill", alpha=.5) +
    geom_text(position="fill", vjust="inward") +
    scale_fill_discrete("Protein has") +
    theme_minimal() +
    scale_x_discrete("") + scale_y_continuous("")

dev.off()

