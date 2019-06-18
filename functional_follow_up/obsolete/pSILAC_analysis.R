##########################################
##
## Pool Silac representation.
##
## For the controls (puf3, puf2, etc.):
##   -line plot with the H/L ratios at 0, 30, 60 and 120 min.
##   -heatmap with the z-scaled H/L ratios at 0, 30, 60 and 120 min.
##   -boxplot FC:2h-0h ~ target/non-target
## 
## For all:
##   -heatmap with genes~ko
##   -boxplot with puf3's targets/non-targets
##   -for the rest of the ko's: boxplot with the targets/non-targets taken from the heatmap.csv analysis (constructs vs. binders)
##
##########################################
options(stringsAsFactors=FALSE)
library(ggplot2)
library(reshape2)
library(rMQanalysis)
library(ggrepel)
library(WriteXLS)
library(gplots)

FC=1
FC2=2  # use to name genes in the lines plot and not overcrowd it

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/")

##
## Analysis of controls
##
pdf("pSILAC_controls.pdf", width=8, height=8)
lapply(c("proteinGroups_adh3_ssd1_pho92.txt",
         "proteinGroups_puf3_puf2_gis2_rad18.txt"), function(f) {

    ##
    ## read input files
    ##
    # read MQ protein groups file
    pg <- read.delim(f)

    # filter out contaminants and filter for unique+razor
    pg_flt   <- filterWholeDataset(pg)
    pg_ident <- filterIdentifiedProteins(pg_flt)

    # add gene name
    ensembl_dictionary <- read.delim("ensembl_dic.txt", sep="\t")
    pg_ident$geneName <- unlist(lapply(strsplit(as.character(pg_ident$Protein.IDs), ";"), function(x) {
      rows <- match(x, ensembl_dictionary$Ensembl.Protein.ID, nomatch=0)
      paste(unique(ensembl_dictionary$Associated.Gene.Name[rows]), collapse=";")
    }))

    pg_ident$geneName <- ifelse(as.character(pg_ident$geneName) == "", 
                                sub("([^;]+);.*", "\\1", pg_ident$Protein.IDs), 
                                sub(";$", "", as.character(pg_ident$geneName)))

    # read the number of different conditions of the experiment 
    conditions <- colnames(pg_ident)[grepl("^Peptides\\..+min_.*", colnames(pg_ident))]
    conditions <- unique(gsub(".+_(.+)$", "\\1", conditions))

    ##
    ##   -line plot with the H/L ratios at 0, 30, 60 and 120 min.
    ##
    lapply(conditions, function(cond) {

      ratio.cols <- paste0(c("Ratio.H.M.normalized.0min_", "Ratio.H.M.normalized.30min_", "Ratio.H.M.normalized.60min_", "Ratio.H.M.normalized.120min_"), cond)

      # keep columns needed for each KO strain of the timecourse
      x <- pg_ident[, c("Protein.IDs", ratio.cols, "geneName")]
      df <- list(proteinID=x$Protein.IDs, 
                 geneName=x$geneName, 
                 quant=log2(as.matrix(x[, grepl("^Ratio", colnames(x))])))
      df$quant[is.na(df$quant)] <- 0
      rownames(df$quant) <- NULL
        
      # organize dataframe for ggplot input
      d <- as.data.frame(df$quant)
      d$gene <- df$geneName
      d <- melt(d)
      d$variable <- factor(d$variable, level=ratio.cols)
        
      # enriched
      genes.enriched <- tapply(d$value, d$gene, function(x) any(abs(x) > FC2))
      d$enriched <- genes.enriched[match(d$gene, names(genes.enriched))]
      d$variable <- ifelse(grepl("^Ratio.H.M.normalized.0min_"  , d$variable), 0,
                    ifelse(grepl("^Ratio.H.M.normalized.30min_" , d$variable), .5,
                    ifelse(grepl("^Ratio.H.M.normalized.60min_" , d$variable), 1,
                    ifelse(grepl("^Ratio.H.M.normalized.120min_", d$variable), 2, -1))))
        
      # plot
      p <- ggplot(d, aes(x=variable, y=value, group=gene, color=enriched)) +
        geom_line() +
        ggtitle(paste0(cond, " KO vs wt")) +
        geom_text_repel(data=subset(d, enriched & (variable == 2)), aes(label=gene), col="black") +
        scale_color_manual(values=c("#00000010", "#FF0000AA")) +
        geom_hline(yintercept = c(-FC, FC), lty=2) +
        geom_hline(yintercept = 0, lty=1) +
        ylim(min(d$value), max(d$value)) +
        theme_minimal() +
        guides(color=FALSE)
        
      print(p)

      invisible(cond)
    })

    ##
    ##   -heatmap with the z-scaled H/L ratios at 0, 30, 60 and 120 min. Only putative targets (>FC)
    ##
    lapply(conditions, function(cond) {

      ratio.cols <- paste0(c("Ratio.H.M.normalized.0min_", "Ratio.H.M.normalized.30min_", "Ratio.H.M.normalized.60min_", "Ratio.H.M.normalized.120min_"), cond)

      # keep columns needed for each KO strain of the timecourse
      x <- pg_ident[, c("Protein.IDs", ratio.cols, "geneName")]
      df <- list(proteinID=x$Protein.IDs, 
                 geneName=x$geneName, 
                 quant=log2(as.matrix(x[, grepl("^Ratio", colnames(x))])))
      df$quant[is.na(df$quant)] <- 0
      rownames(df$quant) <- NULL

      # heatmap
      out <- apply(df$quant, 1, function(x) all(abs(x) < FC))
      heatmap.2(df$quant[!out, ], trace="none", col=redgreen, scale="none", dendrogram="none", Colv=FALSE, main=cond,
                labRow=df$geneName[!out], labCol=c("0min", "30min", "60min", "120min"))

      invisible(cond)
    })
        
    ##
    ##   -boxplot FC:2h-0h ~ target/non-target. Only proteins with at least 1 ratio >0
    ##
    lapply(conditions, function(cond) {

      ratio.cols <- paste0(c("Ratio.H.M.normalized.0min_", "Ratio.H.M.normalized.30min_", "Ratio.H.M.normalized.60min_", "Ratio.H.M.normalized.120min_"), cond)

      # keep columns needed for each KO strain of the timecourse
      x <- pg_ident[, c("Protein.IDs", ratio.cols, "geneName")]
      df <- list(proteinID=x$Protein.IDs, 
                 geneName=x$geneName, 
                 quant=log2(as.matrix(x[, grepl("^Ratio", colnames(x))])))
      df$quant[is.na(df$quant)] <- 0
      rownames(df$quant) <- NULL

      # read targets
      targets <- tryCatch(read.delim(paste0(cond, "_target_mRNAs.txt"))[, 1],
                          error=function(e) { c("") })
      df$target <- sapply(df$proteinID, function(x) any(unlist(strsplit(x, ";")) %in% targets))

      # calculate ratios and prepare the dataset for ggplot2
      df$ratio <- df$quant[, paste0("Ratio.H.M.normalized.120min_", cond)] - df$quant[, paste0("Ratio.H.M.normalized.0min_", cond)] 
      out <- apply(df$quant, 1, function(x) all(x == 0))
      d <- data.frame(ratio=df$ratio[!out],
                      target=df$target[!out])

      # boxplot
      p <- ggplot(d, aes(x=target, y=ratio, color=target)) +
        geom_boxplot() +
        ggtitle(paste(cond, "KO protein abundance ratio at 120min - 0min")) +
        scale_x_discrete("", labels=c("Not Target", "Target")) +
        scale_y_continuous("120 min - 0 min ratio") +
        theme_minimal() +
        guides(color=FALSE)
        
      print(p)

      invisible(cond)
    })
})

dev.off()
  
##
## Analysis of the rest
##
pdf("pSILAC.pdf", width=8, height=8)

# read MQ protein groups file
pg <- read.delim("proteinGroups.txt")

# filter out contaminants and filter for unique+razor
pg_flt   <- filterWholeDataset(pg)
pg_ident <- filterIdentifiedProteins(pg_flt)

# add gene name
ensembl_dictionary <- read.delim("ensembl_dic.txt", sep="\t")
pg_ident$geneName <- unlist(lapply(strsplit(as.character(pg_ident$Protein.IDs), ";"), function(x) {
  rows <- match(x, ensembl_dictionary$Ensembl.Protein.ID, nomatch=0)
  paste(unique(ensembl_dictionary$Associated.Gene.Name[rows]), collapse=";")
}))

pg_ident$geneName <- ifelse(as.character(pg_ident$geneName) == "", 
                            sub("([^;]+);.*", "\\1", pg_ident$Protein.IDs), 
                            sub(";$", "", as.character(pg_ident$geneName)))

# keep columns needed for each KO strain of the timecourse
x <- pg_ident[, grep("(Protein\\.IDs|Ratio\\.H\\.M\\.normalized\\.0min_|Ratio\\.H\\.M\\.normalized\\.120min_|geneName)", colnames(pg_ident))]
df <- list(proteinID=x$Protein.IDs,
           geneName=x$geneName,
           quant=log2(as.matrix(x[, grepl("^Ratio", colnames(x))])))
df$quant[is.na(df$quant)] <- 0
rownames(df$quant) <- NULL

# read the number of different conditions of the experiment and calculate the 120-0 min ratio
conditions <- colnames(pg_ident)[grepl("^Peptides\\..+min_.*", colnames(pg_ident))]
conditions <- unique(gsub(".+_(.+)$", "\\1", conditions))
ratios <- sapply(conditions, function(cond) df$quant[, paste0("Ratio.H.M.normalized.120min_", cond)] - df$quant[, paste0("Ratio.H.M.normalized.0min_", cond)])
rownames(ratios) <- pg_ident$geneName

# heatmap with the H/L ratios at 0 and 120 min. Only putative targets (>FC)
out <- apply(ratios, 1, function(x) all(abs(x) < FC))
d <- ifelse(ratios[!out, ] > FC, 1, ifelse(ratios[!out, ] < -FC, -1, 0))
heatmap.2(d, trace="none", col=redgreen, scale="none", dendrogram="none", Rowv=TRUE, Colv=FALSE, key=FALSE)
write.csv(ratios[!out,], file="pSILAC.csv")

##
## boxplot with puf3 targets
##
# read targets
targets <- tryCatch(read.delim("puf3_target_mRNAs.txt")[, 1])

# remove ratios which the quantification at 0min and 120min are 0
out <- sapply(conditions, function(cond) apply(df$quant[, grepl(cond, colnames(df$quant))], 1, function(x) all(x == 0)))

# boxplot: here we constrain the genes to all that were measures + all targets (regardless they were measured)
is.target   <- sapply(df$proteinID, function(x) any(unlist(strsplit(x, ";")) %in% targets))
is.measured <- !out[, "puf3"]
d <- data.frame(ratio=ratios[is.measured, "puf3"],
                target=is.target[is.measured])

p <- ggplot(d, aes(x=target, y=ratio, color=target)) +
    geom_boxplot() +
    ggtitle("PUF3 KO protein abundance ratio at 120min - 0min") +
    scale_x_discrete("", labels=c("Not Target", paste0("Target (", sum(d$target), ")"))) +
    scale_y_continuous("120 min - 0 min ratio") +
    theme_minimal() +
    guides(color=FALSE)

print(p)

##
##   -for the rest of the ko's: boxplot with the targets/non-targets taken from the heatmap.csv analysis (constructs vs. binders)
##
binders <- read.csv("../analysis/heatmap_1+1_threshold.csv")    # our list of binders <-- the constructs are our putative targets
constructs <- read.csv("../constructs/constructs.csv")  # translate the construct name from P1A1 to gene name

ratios2 <- ratios[, colnames(ratios) != "puf3"]   # remove puf3 (positive control, we already know its binders)
out2    <- out   [, colnames(out   ) != "puf3"]

sapply(1:ncol(ratios2), function(i) {

    binder <- toupper(colnames(ratios2)[i])

    if(any(binders$name == binder)) {
        targets <- binders[binders$name == binder, c(-1, -2)] > 1 # which constructs enrich this binder?
        colnames(targets) <- constructs$gene[match(colnames(targets), constructs$name)]  # construct names to gene names

        # boxplot: here we constrain the genes to all that were measured + all targets (regardless they were measured)
        is.target   <- sapply(df$proteinID, function(x) any(unlist(strsplit(x, ";")) %in% colnames(targets)[targets]))
        is.measured <- !out2[, i]
        d <- data.frame(ratio=ratios2[is.measured, i],
                        target=is.target[is.measured])
        p <- ggplot(d, aes(x=target, y=ratio, color=target)) +
            geom_boxplot() +
            ggtitle(paste(colnames(ratios2)[i], "KO protein abundance ratio at 120min - 0min")) +
            scale_x_discrete("", labels=c("Not Target", paste0("Target (", sum(d$target), ")"))) +
            scale_y_continuous("120 min - 0 min ratio") +
            theme_minimal() +
            guides(color=FALSE)

        print(p)
    }

    invisible(colnames(ratios2)[i])
})

dev.off()
