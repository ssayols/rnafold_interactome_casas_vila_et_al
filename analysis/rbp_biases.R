#####################################
##
## On RBPs list check for bias towards: 
##   * larger → db/protein_properties.SGD.tab
##   * more hydrophobic → db/protein_properties.SGD.tab
##   * more basic pH protein identification → db/protein_properties.SGD.tab
##   * Disordered regions → get domain/disorder data from http://pfam.xfam.org,
##     and calculate fraction of disorder for each protein. Plot per RBP vs. bg
##     fraction of proteins vs. fraction of disorder.Also, repeat separating
##     proteins with known RBD (pfam) versus unknown.
##
## See description in of the features in https://sites.google.com/view/yeastgenome-help/function-help/protein-information
##
#####################################
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(rMQanalysis)

set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/rbp_biases.pdf", onefile=TRUE)

##
## 1-read input data
##
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv")    # our list of binders
bg      <- filterIdentifiedProteins(filterWholeDataset(read.delim("./analysis/MQ/proteinGroups_reference_lysate.txt")))

# deduplicate binders
binders <- by(binders, 1:nrow(binders), function(x) {
    data.frame(X=unlist(strsplit(x$X, ";")),
               x[,!grepl("^(X|name)", colnames(x))],    # gene name is no longer valid for XXX;YYY rows
               row.names=NULL)
})
binders <- do.call(rbind, binders)

# deduplicate bg
bg <- do.call(rbind, by(bg, 1:nrow(bg), function(x) {
    data.frame(X=unlist(strsplit(x$Majority.protein.IDs, ";")),
               x[,colnames(x) != "Majority.protein.IDs"],
               row.names=NULL)
}))

# the second background: all yeast proteins (measured or not)
bg2 <- getBM(attributes=c("ensembl_gene_id", "peptide"),
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
bg2 <- unique(bg2, MARGIN=1)
bg2 <- bg2[bg2$peptide != "bg2uence unavailable", ]
bg2 <- bg2[sapply(bg2$peptide, nchar) > 70, ]  # remove proteins shorter than 70aa (~250)
colnames(bg2)[1] <- "X"

# read protein properties and calculate fraction of aa. instead of absolute counts
protein_properties <- read.delim("./analysis/db/protein_properties.SGD.tab")
protein_properties <- protein_properties[protein_properties$ORF %in% unique(c(bg$X, bg2$X, binders$X)), ]

# some transformations
protein_properties$Mw <- log2(protein_properties$Mw)
protein_properties$Protein.Length <- log2(protein_properties$Protein.Length)
protein_properties$ASSUMING.ALL.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES <- log2(protein_properties$ASSUMING.ALL.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES)
protein_properties$ASSUMING.NO.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES <- log2(protein_properties$ASSUMING.NO.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES)
aa <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
protein_properties[, aa] <- t(apply(protein_properties[, aa], 1, function(x) x / sum(x)))

# biases of interest
cols <- c("Mw", "PI", "Protein.Length", "GRAVY.Score", "Aromaticity.Score", "CAI",
          "Codon.Bias", "FOP.Score", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly",
          "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser",
          "Thr", "Val", "Trp", "Tyr", "CARBON", "HYDROGEN", "NITROGEN", "OXYGEN",
          "SULPHUR", "INSTABILITY.INDEX..II.", "ASSUMING.ALL.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES",
          "ASSUMING.NO.CYS.RESIDUES.APPEAR.AS.HALF.CYSTINES", "ALIPHATIC.INDEX")

##
## 2-density of proteins vs. each confounder in binders and bg 
##
Map(function(x, col, orf=protein_properties$ORF) {
  x <- data.frame(x=x)
  i <- !is.infinite(x$x) & !is.na(x$x)
  orf <- orf[i]
  x <- x[i, , drop=FALSE]
  x$class1 <- ifelse(orf %in% binders$X, "identified",
              ifelse(orf %in% bg$X, "measured", "predicted_only"))
  x$class2 <- ifelse(orf %in% binders$X, "identified", "predicted")
  pval_id_me <- with(x, format.pval(t.test(x[class1 == "identified"], x[class1 == "measured" ])$p.value))
  pval_me_pr <- with(x, format.pval(t.test(x[class1 == "measured"  ], x[class2 == "predicted"])$p.value))

  ggplot() +
    geom_density(data=subset(x, class1 %in% c("identified", "measured")), mapping=aes(x=x, color=class1, fill=class1), alpha=.5) +
    geom_density(data=subset(x, class2 %in% c("predicted"))             , mapping=aes(x=x, color=class2, fill=class2), alpha=.5) +
    labs(title=col,
         subtitle=paste0("identified vs. measured (p=", pval_id_me, ")",
                         " measured vs. predicted (p=", pval_me_pr, ")"),
         caption="check description in https://sites.google.com/view/yeastgenome-help/function-help/protein-information",
         x=col) +
    scale_color_discrete("") +
    scale_fill_discrete("") +
    theme_bw()
}, protein_properties[, cols], cols)

dev.off()
