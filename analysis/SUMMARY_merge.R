#########################
##
## merge data from several analyses tables into the SUMMARY excel shit
##
#########################
library(readxl)
library(WriteXLS)
library(biomaRt)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")

ALL <- read_excel("./analysis/SUMMARY (Nuria).xlsx", sheet="ALL")

##
## number of constructs the binder binds to
##
hm <- read.csv("./analysis/heatmap_1+1_threshold.csv")

# add ensembl_gene_id
ALL$yORF <- hm$X[match(ALL$protein, hm$name)]

# count in how many constructs each protein binds to (R1 criteria)
x <- apply(hm[, c(-1, -2)], 1, function(x) sum(x > 1))
ALL$"binds to constructs" <- x[match(ALL$yORF, hm$X)]

# split rows if there were multiple genes per row
ALL <- apply(ALL, 1, function(x) {
    data.frame("protein"=x["protein"],
               "whole name"=x["whole name"],
               "biological process/function"=x["biological process/function"],
               "RNA binding (protein type)"=x["RNA binding (protein type)"],
               "RNA type"=x["RNA type"],
               "localization"=x["localization"],
               "description"=x["description"],
               "phenotye"=x["phenotye"],
               "yORF"=unlist(strsplit(x["yORF"], ";")),
               "binds to constructs"=x["binds to constructs"])
})
ALL <- do.call(rbind, ALL)

##
## localization data from Weissman
##
loc <- read.delim("./analysis/db/localization.txt", comment.char="#")
ALL$compartment <- sapply(ALL$yORF, function(x) {
    paste(sapply(unlist(strsplit(x, ";")), function(x) {
        x <- sapply(c("nucleus", "nucleolus", "cytoplasm", "mitochondrion"), function(y) {
            if(loc[loc$yORF == x, y]) y else ""
        })
        paste(x[x != ""], collapse=", ")
    }), collapse="; ")
})

##
## GO data for binders
##
rbd <- read.csv("./analysis/domains_binders.csv")
ALL$go_rna_bind <- rbd$go_rna_bind[match(ALL$yORF, rbd$X)]
ALL$go_rna_rel  <- rbd$go_rna_rel[match(ALL$yORF, rbd$X)]

##
## RBP data for binders, taken from the Tuschel's census
##
pfam     <- read.delim("./analysis/db/Pfam-A.clans.tsv", head=F)    # translation between ensembl pfam_id and Tuschel's pfam name
rbd.pfam <- read.delim("./analysis/db/PfamA_tuschel.csv")
rbd.go    <- read.delim("./analysis/db/GO_tuschel.csv")
 
## COPIED FROM DOMAINS.R
annotate.domains <- function(x, dataset, key) {

    # get gene's pfam domains and GO terms from biomart
    ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "pfam", "name_1006"),
                      filters=c("ensembl_gene_id"), values=x[, key],
                      mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset=dataset))

    # check which x are RNA related from Tuschel's pfamA list
    x$pfam <- sapply(x[, key], function(x) {
        pfam_id   <- ann$pfam[ann$ensembl_gene_id == x]
        pfam_name <- pfam$V4[match(pfam_id, pfam$V1, nomatch=0)]
        paste(unique(pfam_name[pfam_name %in% rbd.pfam$Pfam.RNA.binding.domains]), collapse=", ")
    })

    # check which x have the "RNA binding" GO term
    x$go_rna_bind <- sapply(x[, key], function(x) any(grepl("\\bRNA\\b binding", ann$name_1006[ann$ensembl_gene_id == x])))

    # check which x are RNA related from Tuschel's GO list
    x$go_rna_rel <- sapply(x[, key], function(x)
        paste(unique(ann$name_1006[ann$name_1006[ann$ensembl_gene_id == x] %in% rbd.go$selected.GO.terms.for.RN..related.processes]),
              collapse=", "))

    x
}
ALL <- annotate.domains(ALL, "scerevisiae_gene_ensembl", "yORF")

##
## which human ortholog + is the ortholog rbp?
##
## COPIED FROM DOMAINS.R
# get human orthologs of the binders with non RNA related domains
orth <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene",
                           "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_orthology_confidence"),
              filters=c("ensembl_gene_id"), values=ALL$yORF,
              mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))

# annotate the human domains
colnames(orth)[grepl("hsapiens_homolog_ensembl_gene", colnames(orth))] <- "yORF"
orth <- annotate.domains(orth, "hsapiens_gene_ensembl", "yORF")

# collapse rows by saccer ensembl gene id
orth <- by(orth, orth$ensembl_gene_id, function(x) {
    data.frame(ensembl_gene_id=paste(unique(x$ensembl_gene_id), collapse=";"),
               external_gene_name=paste(unique(x$external_gene_name), collapse=";"),
               hsapiens_homolog_ensembl_gene=paste(unique(x$yORF), collapse=";"),
               hsapiens_homolog_associated_gene_name=paste(unique(x$hsapiens_homolog_associated_gene_name), collapse=";"),
               pfam=x$pfam,
               go_rna_bind=any(x$go_rna_bind),
               go_rna_rel=any(x$go_rna_rel != ""))
})
orth <- do.call(rbind, orth)
orth <- unique(orth, MARGIN=1)

ALL$hsapiens_homolog             <- orth$hsapiens_homolog_associated_gene_name[match(ALL$yORF, orth$ensembl_gene_id)]
ALL$hsapiens_homolog_pfam        <- orth$pfam[match(ALL$yORF, orth$ensembl_gene_id)]
ALL$hsapiens_homolog_go_rna_bind <- orth$go_rna_bind[match(ALL$yORF, orth$ensembl_gene_id)]
ALL$hsapiens_homolog_go_rna_rel  <- orth$go_rna_rel[match(ALL$yORF, orth$ensembl_gene_id)]

##
## is Mitchell's table?
##
mitchell <- read.delim("./analysis/db/yeastRBP_Mitchell.txt", comment="#")  # former Parker
ALL$Mitchell <- ALL$yORF %in% mitchell$Systematic.Name

##
## in Matia_Gonzalez's table?
##
matia_gonzalez <- read.delim("./analysis/db/yeastRBP_Matia_Gonzalez.txt", comment="#")  # former Gerber
ALL$Matia_Gonzalez <- ALL$yORF %in% matia_gonzalez$ORF

##
## in Merve's table? take ProteinID col
##
merve <- read.delim("./analysis/db/Merve.txt")
ALL$Merve <- ALL$yORF %in% unlist(strsplit(merve$Protein.IDs, ";"))

##
## in haploid KO table? ORF.name col
##
haploid <- read_excel("./analysis/db/KO_haploid.xls", sheet="mat_alpha_obs")
ALL$KOhaploid <-  ALL$yORF %in% haploid$"ORF name"

##
## in overexpression table? Sys. Name col
##
oe <- read_excel("./analysis/db/overexpression.xls", sheet="Sheet1")
ALL$overexpression <-  ALL$yORF %in% oe$"Sys. Name"

##
## in which filtering criterias we identify the binder? r1, rpuf3 or r1+1
##
hm.1    <- read.csv("./analysis/heatmap_1+1_threshold.csv")
hm.puf3 <- read.csv("./analysis/heatmap_puf3_threshold.csv")

ALL$"r1 filter"    <- ALL$yORF %in% unlist(strsplit(hm$X, ";"))
ALL$"r1+1 filter"  <- ALL$yORF %in% unlist(strsplit(hm.1$X, ";"))
ALL$"rpuf3 filter" <- ALL$yORF %in% unlist(strsplit(hm.puf3$X, ";"))

##
## save results back
##
l <- list(ALL)
WriteXLS("l", ExcelFileName="./analysis/SUMMARY.xls", SheetNames=c("Sheet1"))
