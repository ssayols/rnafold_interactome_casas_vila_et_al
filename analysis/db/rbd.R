########################
##
## Annotate a list of yeast proteins with rna-binding domains (from Tuschel's census)
##
########################
library(rMQanalysis)
library(biomaRt)
library(WriteXLS)

bg1 <- filterIdentifiedProteins(filterWholeDataset(read.delim("./analysis/MQ/proteinGroups_reference_lysate.txt")))
bg2 <- getBM(attributes=c("ensembl_gene_id"), mart=useMart("ensembl", dataset="scerevisiae_gene_ensembl"))
colnames(bg) <- "X"

pfam     <- read.delim("./analysis/db/Pfam-A.clans.tsv", head=F)    # translation between ensembl pfam_id and Tuschel's pfam name
rbd.pfam <- read.delim("./analysis/db/PfamA_tuschel.csv")

out <- lapply(list(llisat=bg1, transcriptome=bg2), function(bg) {
    ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "pfam", "name_1006"),
                 filters=c("ensembl_gene_id"), values=unlist(strsplit(bg$X, ";")),
                 mart=useMart("ensembl",dataset="scerevisiae_gene_ensembl"))

    bg$rbp <- sapply(strsplit(bg$X, ";"), function(x) {
        pfam_id   <- ann$pfam[ann$ensembl_gene_id %in% x]
        pfam_name <- pfam$V4[match(pfam_id, pfam$V1, nomatch=0)]
        any(pfam_name %in% rbd.pfam$Pfam.RNA.binding.domains)
    })

    bg$gene <- sapply(strsplit(bg$X, ";"), function(x) {
        paste(unique(ann$external_gene_name[ann$ensembl_gene_id %in% x]), collapse=";")
    })
})

WriteXLS("out",ExcelFileName="./analysis/db/rbd.xls")
