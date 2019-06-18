options(stringsAsFactors=FALSE)
library(ggplot2)
library(reshape2)
library(rMQanalysis)
library(ggrepel)
library(WriteXLS)

##
## read protein groups file
##
# convert ensembl yeast IDs to ensembl Associated gene names (using biomart 'dicctionary')
ensembl_dictionary <- read.delim("ensembl_dic.txt", sep="\t")
pg <- read.delim("proteinGroups_puf3_puf2_gis2_rad18.txt")
pg_ensembl <- strsplit(as.character(pg$Protein.IDs), ";")

pg$Gene.Name <- unlist(lapply(pg_ensembl, function(x) {
  rows <- match(x, ensembl_dictionary$Ensembl.Protein.ID, nomatch=0)
  paste(unique(ensembl_dictionary$Associated.Gene.Name[rows]), collapse=";")
}))

# filter out contaminants and filter for unique+razor
pg_flt   <- filterWholeDataset(pg)
pg_ident <- filterIdentifiedProteins(pg_flt)

# create a $my_label column
pg_ident$my_label <- ifelse(as.character(pg_ident$Gene.Name) == "", 
                            sub("([^;]+);.*", "\\1", pg_ident$Protein.IDs), 
                            sub(";$', '", as.character(pg_ident$Gene.Name)))

# read the number of different conditions of the experiment 
conditions <- colnames(pg_ident)[grepl("^Peptides\\..+min_.*", colnames(pg_ident))]
conditions <- unique(gsub(".+_(.+)$", "\\1", conditions))

##
## plot silac ratios along time
##
pdf("pSILAC_ggplot2.pdf", width=8, height=8)
results_ggplot <- lapply(conditions, function(cond) {

  #keep columns needed for each KO strain of the timecourse
  x <- pg_ident[, c("Protein.IDs", paste0(c("Ratio.H.M.0min_", "Ratio.H.M.30min_", "Ratio.H.M.60min_", "Ratio.H.M.120min_"), cond), "my_label")]
  df <- list(proteinID=x$Protein.IDs, 
             geneName=x$my_label, 
             quant=log2(as.matrix(x[, grepl("^Ratio", colnames(x))])))
  df$quant[is.na(df$quant)] <- 0
  rownames(df$quant) <- NULL
    
  # organize dataframe for ggplot input
  d <- as.data.frame(df$quant)
  d$gene <- df$geneName
  d <- melt(d)
  d$variable <- factor(d$variable, level=paste0(c("Ratio.H.M.0min_", "Ratio.H.M.30min_", "Ratio.H.M.60min_", "Ratio.H.M.120min_"), cond))
    
  # enriched
  genes.enriched <- tapply(d$value, d$gene, function(x) any(abs(x) > 2))
  d$enriched <- genes.enriched[match(d$gene, names(genes.enriched))]
  d$variable <- ifelse(d$variable == "Ratio.H.M.0min_adh3"  , 0,
                ifelse(d$variable == "Ratio.H.M.30min_adh3" , .5,
                ifelse(d$variable == "Ratio.H.M.60min_adh3" , 1,
                ifelse(d$variable == "Ratio.H.M.120min_adh3", 2, -1))))
    
  # plot
  p <- ggplot(d, aes(x=variable, y=value, group=gene, color=enriched)) +
    geom_line() +
    ggtitle(paste0(cond, " KO vs wt")) +
    geom_text_repel(data=subset(d, enriched & variable == paste0("Ratio.H.M.120min_", cond)), aes(label=gene), col="black") +
    scale_color_manual(values=c("#00000010", "#FF0000AA")) +
    geom_hline(yintercept = c(-1, 1), lty=2) +
    geom_hline(yintercept = 0, lty=1) +
    ylim(min(d$value), max(d$value)) +
    theme_minimal() +
    guides(color=FALSE)
    
  print(p)
})
  
dev.off()
