##########################################
##
## Pool Silac representation.
##
##  -heatmap with genes~ko
##  -boxplot with the targets/non-targets taken from the heatmap.csv analysis (constructs vs. binders)
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

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/")
pdf("pSILAC_with_names.pdf", width=4, height=4)

##
## Read pSILAC data, and get the 0-120 H/M ratios
##
# read known targets from controls
known_targets <- lapply(list(GIS2="../functional_follow_up/gis2_target_mRNAs.txt",
                             PHO92="../functional_follow_up/pho92_target_mRNAs.txt",
                             PUF2="../functional_follow_up/puf2_target_mRNAs.txt",
                             PUF3="../functional_follow_up/puf3_target_mRNAs.txt",
                             SSD1="../functional_follow_up/ssd1_target_mRNAs.txt"),
                        function(f) read.delim(f)[, 1])

# read MQ protein groups file and filter out contaminants and filter for unique+razor
pg <- read.delim("proteinGroups_pSILAC.txt")
pg_flt   <- filterWholeDataset(pg)
pg_ident <- filterIdentifiedProteins(pg_flt)

# add gene name
ensembl_dictionary <- read.delim("ensembl_dic.txt", sep="\t")
pg_ident$geneName <- sapply(strsplit(pg_ident$Protein.IDs, ";"), function(x) {
  rows <- match(x, ensembl_dictionary$Ensembl.Gene.ID, nomatch=0)
  paste(unique(ensembl_dictionary$Associated.Gene.Name[rows]), collapse=";")
})

pg_ident$geneName <- ifelse(as.character(pg_ident$geneName) == "", 
                            sub("([^;]+);.*", "\\1", pg_ident$Protein.IDs), 
                            sub(";$", "", as.character(pg_ident$geneName)))

# read the number of different conditions of the experiment
conditions <- colnames(pg_ident)[grepl("^Peptides\\.0min_.*", colnames(pg_ident))]
conditions <- sort(unique(gsub(".+?_(.+)$", "\\1", conditions)))

# keep columns needed for each KO strain of the timecourse
x <- pg_ident[, grep("(Protein\\.IDs|geneName|Ratio\\.H\\.M\\.normalized\\.0min_|Ratio\\.H\\.M\\.normalized\\.120min_)", colnames(pg_ident))]
df <- list(proteinID=pg_ident$Protein.IDs,
           geneName=pg_ident$geneName,
           ratio0  =log2(as.matrix(pg_ident[, paste0("Ratio.H.M.normalized.0min_"  , conditions)])),
           ratio120=log2(as.matrix(pg_ident[, paste0("Ratio.H.M.normalized.120min_", conditions)])))
df$ratio0  [is.na(df$ratio0  )] <- 0
df$ratio120[is.na(df$ratio120)] <- 0
rownames(df$ratio0  ) <- NULL
rownames(df$ratio120) <- NULL
colnames(df$ratio0  ) <- gsub("^Ratio.H.M.normalized.0min_"  , "", colnames(df$ratio0  ))
colnames(df$ratio120) <- gsub("^Ratio.H.M.normalized.120min_", "", colnames(df$ratio120))

is.replicated <- sapply(conditions, function(x) regexpr(".+_(\\d)$", x)[[1]] > 0)
is.measured   <- apply(df$ratio120, 2, function(x) abs(x) > 0)

##
## heatmap with the ratios. Only putative targets (>FC) in at least 1 condition
##
out <- apply(df$ratio120, 1, function(x) all(abs(x) < FC))
d   <- df$ratio120[!out, ]
heatmap.2(d, trace="none", col=function(n) colorpanel(n, "red", "white", "green"), scale="none",
          dendrogram="none", Rowv=TRUE, Colv=TRUE, key=TRUE, colsep=1:ncol(d), sepcolor="black")
d   <- ifelse(df$ratio120[!out, ] >  FC,  1,
       ifelse(df$ratio120[!out, ] < -FC, -1, 0))
heatmap.2(d, trace="none", col=function(n) colorpanel(n, "red", "white", "green"), scale="none",
          dendrogram="none", Rowv=TRUE, Colv=TRUE, key=TRUE, colsep=1:ncol(d), sepcolor="black")

##
##   boxplot with the targets/non-targets taken from the heatmap.csv analysis (constructs vs. binders)
##
binders    <- read.csv("../analysis/heatmap_1+1_threshold.csv")    # our list of binders <-- the constructs are our putative targets
constructs <- read.csv("../constructs/constructs.csv")  # translate the construct name from P1A1 to gene name

res <- lapply(1:ncol(df$ratio120), function(i) {

    binder <- toupper(colnames(df$ratio120)[i])
    if(is.replicated[i]) binder <- gsub("_\\d+$", "", binder)

    # if it's in our binders table or is a control (puf2, puf3, etc.), process
    if(any(j <- grepl(binder, binders$name, fixed=TRUE)) | binder %in% names(known_targets)) {

        targets.of.this.binder <-
            if(binder %in% names(known_targets)) {  # 'i' is a control
                known_targets[[binder]]
            } else {  # 'i' is in our binders table
                targets <- abs(binders[j, c(-1, -2)]) > 1  # which constructs enrich this binder?
                constructs$gene[match(colnames(targets)[targets], constructs$name)]  # construct names to gene IDs
            }

        # boxplot: here we constrain the genes to all that were measured
        is.target <- sapply(strsplit(df$proteinID, ";"), function(x) any(x %in% targets.of.this.binder))
        d <- data.frame(proteinID=df$proteinID[is.measured[, i]],
                        geneName =df$geneName [is.measured[, i]],
                        ratio0   =df$ratio0   [is.measured[, i], i],
                        ratio120 =df$ratio120 [is.measured[, i], i],
                        target   =is.target[is.measured[, i]],
                        binder   =colnames(df$ratio120)[i])
        p <- ggplot(d, aes(x=target, y=ratio120, color=target)) +
                geom_boxplot() +
                geom_text_repel(data=subset(d, abs(d$ratio120) > FC), mapping=aes(x=target, y=ratio120, label=geneName)) +
                ggtitle(paste(colnames(df$ratio120)[i], "KO protein abundance ratio at 120min")) +
                scale_x_discrete("", labels=c(paste0("Not Target (", sum(!d$target), ")"), paste0("Target (", sum(d$target), ")"))) +
                scale_y_continuous("120 min H/M ratio") +
                theme_minimal() +
                guides(color=FALSE)
        print(p)

        d[d$target | abs(d$ratio120) > FC, ]
    } else {
        data.frame(proteinID="", geneName="", ratio0=0, ratio120=0, target=F, binder=colnames(df$ratio120)[i])
    }
})

write.csv(do.call(rbind, res), file="pSILAC.csv", row.names=F)
dev.off()
