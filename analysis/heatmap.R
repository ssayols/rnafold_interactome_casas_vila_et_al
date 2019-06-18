#############################################
##
## analysis of the proteins bound to the constructs
##
#############################################
options(stringsAsFactors=F)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(rMQanalysis)

CONTROLS <- c(LSG1="YGL099W", SUI3="YPL237W", GCD11="YER025W", PUF3="YLL013C")
PG <- "./analysis/MQ/proteinGroups.txt"
CONSTRUCTS <- "./constructs/constructs.csv"
GOSLIM <- "./analysis/db/go_slim_mapping.tab"
OUT <- "./analysis"

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf(paste0(OUT, "/heatmap_1+1_threshold.pdf"), onefile=TRUE, width=14, height=14)

##
## read input data, filter out contaminants and get the ratio columns only
##
pg_filtered <- filterIdentifiedProteins(filterWholeDataset(read.delim(PG)))

# construir dataframe amb els ratios
keep_column <- colnames(pg_filtered)[grepl("Ratio\\.H\\.L\\.normalized\\.P.+", colnames(pg_filtered))]
experiments <- unique(gsub("Ratio\\.H\\.L\\.normalized\\.(P.+)_.+", "\\1", keep_column))
df_ratio <- pg_filtered[, c("Majority.protein.IDs", keep_column)]

# construir dataframe amb les columnes
keep_column <- colnames(pg_filtered)[grepl("(Intensity\\.H\\.P.+_for$|Intensity\\.L\\.P.+_rev$)", colnames(pg_filtered))]
df_intensity <- pg_filtered[, c("Majority.protein.IDs", keep_column)]

##
## some QC: barplot with SE bars of the PUF3 forward and reverse experiments for all constructs
##
f <- log2(as.numeric(df_ratio[df_ratio$Majority.protein.IDs == CONTROLS["PUF3"], grepl("_for$", colnames(df_ratio))]))
r <- log2(as.numeric(df_ratio[df_ratio$Majority.protein.IDs == CONTROLS["PUF3"], grepl("_rev$", colnames(df_ratio))]))
v <- sqrt(f^2 + r^2)
pal <- RColorBrewer::brewer.pal(3, "Set2")
boxplot(data.frame(forward=abs(f), reverse=r, vector=v), col=pal, main="distribution of PUF3 SILAC ratios", bty="n")
#err <- c(sd(f) / sqrt(length(f)), sd(r) / sqrt(length(r)))
sdev <- c(sd(f, na.rm=T), sd(r, na.rm=T), sd(v, na.rm=T))
avg  <- c("|forward|"=mean(f, na.rm=T), reverse=mean(r, na.rm=T), vector=mean(v, na.rm=T))
b <- barplot(abs(avg), col=pal, ylim=c(0, max(avg + sdev)), main="distribution of PUF3 SILAC ratios")
arrows(b, abs(avg) + sdev, b, abs(avg), angle=90, code=1, length=.05)

##
## calculated which proteins are enriched
##
df_ratio_mean <- sapply(experiments, function(e){
    forward <- paste0("Ratio.H.L.normalized.", e, "_for")
    reverse <- paste0("Ratio.H.L.normalized.", e, "_rev")
    x <- log2(df_ratio[, reverse])
    y <- log2(df_ratio[, forward])
    
    # calcular h:
    #   -qualitativament
    #   -nomes al upper-left (binders) i lower-right quadrants (controls)
    #   -log ratios als eixos x i y > RATIO
    x <- ifelse(is.na(x), 0, x)
    y <- ifelse(is.na(y), 0, y)

    # segons: #1 --> el control de referencia (PUF3)
    #     --> #2 --> threshold circular de radi 1 al voltant de l'origen (0,0)
    #         #3 --> threshold lineal a x=1 i y=1
#   min.th  <- 0                                                               #12
    min.th  <- 1                                                               #3
#   th.puf3.x <- x[df_ratio$Majority.protein.IDs == CONTROLS["PUF3"]]          #1
#   th.puf3.y <- y[df_ratio$Majority.protein.IDs == CONTROLS["PUF3"]]          #1
#   THRESHOLD <- sqrt(th.puf3.x^2 + th.puf3.y^2)  # mida vector del control    #1
    THRESHOLD <- 1                                                             #23
    vlen  <- sqrt(x^2 + y^2)    # vector length of the forward and reverse experiments
    ratio <- vlen / THRESHOLD   # enrichment with respect to the threshold
    h <- ifelse(x < -min.th & y >  min.th & vlen >= THRESHOLD,  ratio, # binder side 
         ifelse(x >  min.th & y < -min.th & vlen >= THRESHOLD, -ratio, # control side
                                                                    0))# bg          

    h <- ifelse(!(df_ratio$Majority.protein.IDs %in% CONTROLS) & h < -1, 0, h) # desmarcar els no controls
    h
})
df_ratio_mean <- as.data.frame(df_ratio_mean)
rownames(df_ratio_mean) <- df_ratio[ ,1]
colnames(df_ratio_mean) <- experiments

##
## and plot heatmaps
##
m <- as.matrix(df_ratio_mean[!apply(df_ratio_mean, 1, function(x) all(x == 0)), ])   # agafar nomes els binders/controls
m <- m[order(rowSums(m)), ]
ids.m <- rownames(m)

# traduir nom de les proteines UNIPROT --> Wikigene name
ann <- getBM(attributes=c("external_gene_name", "ensembl_gene_id"), 
             filters="ensembl_gene_id", values=unique(unlist(strsplit(ids.m, ";"))), 
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
names.m <- sapply(ids.m, function(x) {
    ensembl_genes  <- unlist(strsplit(x, ";"))
    external_genes <- ann$external_gene_name[match(ensembl_genes, ann$ensembl_gene_id)]
    if(all(external_genes == "")) x else paste(unique(external_genes), collapse=";")
})

# posar color al constructe segons la posicio
constructs <- read.csv(CONSTRUCTS)
constructs$name <- paste0(constructs$plate, constructs$position)
constructs$pos  <- gsub(".+_(.+)_.$", "\\1", constructs$name1)

fer.heatmap <- function(m, rown, pal, RowSideColors, main="Proteins with higher SILAC ratios than 1") {
  if(pal == "yellows") {
    pal <- c(colorRampPalette(c(rep("#E41A1C", 2), "black"))(50),
             colorRampPalette(c("black", rep("yellow", 2)))(50))
  } else if(pal == "blues") {
    pal <- c(colorRampPalette(c(rep("#E41A1C", 2), "white"))(50),
             colorRampPalette(c("white", rep("#377EB8", 2)))(50))
  } else {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, pal))(100)
  }
  pos <- factor(constructs$pos[match(colnames(m), constructs$name)])
  cc  <- brewer.pal(3, "Set1")[pos]
  o <- order(pos, colSums(m))
  h <- heatmap.2(m[, o], ColSideColors=cc[o], 
#                 colsep=1:ncol(m), rowsep=1:nrow(m), sepwidth=c(.01, .01),
                 RowSideColors=RowSideColors, trace="none", 
                 Rowv=T, Colv=F, col=pal, dendrogram="none", cexRow=0.5, cexCol=0.5, key=T, 
                 labRow=rown, labCol=constructs$gene[match(colnames(m), constructs$name)],
                 density.info="none", main=main)
  #legend("topright", fill=brewer.pal(3, "Pastel1"), legend=levels(pos))
  invisible(h)
}

# totes les proteines
unique.binders <- apply(m, 1, function(x) sum(x > 0) == 1)
h <- fer.heatmap(m=m, rown=names.m, pal="RdBu", RowSideColors=brewer.pal(3, "Set1")[as.factor(unique.binders)])
h <- fer.heatmap(m=m, rown=names.m, pal="yellows", RowSideColors=brewer.pal(3, "Set1")[as.factor(unique.binders)])

out <- as.data.frame(t(h$carpet))
out$name <- names.m[match(rownames(out), names(names.m))]
out <- out[, order(colnames(out))]
write.csv(out, file=paste0(OUT, "/heatmap_1+1_threshold.csv"))

# nomes els unique binders
fer.heatmap(m=m[unique.binders, ], rown=names.m[unique.binders], pal="Blues", RowSideColors=rep("#CCEBC5", sum(unique.binders)))

# nomes els common binders
fer.heatmap(m=m[!unique.binders, ], rown=names.m[!unique.binders], pal="RdBu", RowSideColors=rep("#FBB4AE", sum(!unique.binders)))

# nomes els controls
fer.heatmap(m=m[CONTROLS, ], rown=names.m[!unique.binders], pal="RdBu")

##
## els mateixos 3 heatmaps, pero treient les proteines ribosomals
##
# llegir i filtrar els GO slim de SGD
goslim <- read.delim(GOSLIM, comment.char="#")
goslim.ribo <- goslim[grepl("ribo|rrna", goslim$term, ignore.case=T),]
ribo <- names.m %in% goslim.ribo$gene

# totes les proteines
h <- fer.heatmap(m=m[!ribo,], rown=names.m[!ribo],
                 pal="RdBu", main="Non-ribosomal proteins\nwith higher SILAC ratios than 1",
                 RowSideColors=c("#FBB4AE", "#CCEBC5")[as.factor(unique.binders[!ribo])])

out <- as.data.frame(t(h$carpet))
out$name <- names.m[match(rownames(out), names(names.m))]
out <- out[, order(colnames(out))]
write.csv(out, file=paste0(OUT, "/heatmap_1+1_threshold.noRibo.csv"))

# nomes els unique binders
h <- fer.heatmap(m=m[unique.binders & !ribo,], rown=names.m[unique.binders & !ribo],
                 pal="Blues", main="Non-ribosomal proteins\nwith higher SILAC ratios than 1",
                 RowSideColors=rep("#CCEBC5", sum(unique.binders & !ribo)))

# nomes els common binders
h <- fer.heatmap(m=m[!unique.binders & !ribo,], rown=names.m[!unique.binders & !ribo],
                 pal="RdBu", main="Non-ribosomal proteins\nwith higher SILAC ratios than 1",
                 RowSideColors=rep("#FBB4AE", sum(!unique.binders & !ribo)))

##
## comprovar que no saturem el sistema amb proteines ribosomals
##
f <- function(e, id){
    forward <- paste0("Intensity.H.", e, "_for")
    reverse <- paste0("Intensity.L.", e, "_rev")
    x <- log(df_intensity[df_intensity$Majority.protein.IDs %in% id, reverse])
    y <- log(df_intensity[df_intensity$Majority.protein.IDs %in% id, forward])
    x <- ifelse(is.infinite(x), 0, x)
    y <- ifelse(is.infinite(y), 0, y)
    mapply(function(x, y) mean(c(x, y), na.rm=T), x, y)
}
df_intensities <- sapply(experiments, f, id=ids.m)
rownames(df_intensities) <- df_intensity$Majority.protein.IDs[df_intensity$Majority.protein.IDs %in% ids.m]
df_intensities <- df_intensities[match(ids.m, rownames(df_intensities)), ]
df_random      <- sapply(experiments, f, id=sample(df_intensity$Majority.protein.IDs, 100))

# boxplot dels enriquiments, separant les proteines per quants constructes s'enganxen
group <- apply(m, 1, function(x) sum(x > 0))
barplot(table(group), ylab="Number of proteins in this category", xlab="Number of constructs the protein binds to")
group <- ifelse(group > 1, 2, group)
barplot(table(group), col=RColorBrewer::brewer.pal(length(table(group)), "Set2"),
        ylab="Number of proteins in this category", xlab="Number of constructs the protein binds to")
enrichments <- sapply(sort(unique(group)), function(g) { # mitja de les intensitats als constructes on la proteina esta enriquida
    unlist(lapply(which(group == g), function(i) {
        mean(df_intensities[i, m[i, ] != 0 & !is.na(m[i, ])], na.rm=T)
    }))
})
names(enrichments) <- sort(unique(group))
enrichments$random <- apply(df_random, 1, mean, na.rm=T)
pal <- RColorBrewer::brewer.pal(length(enrichments), "Set2")
boxplot(enrichments, col=pal, bty="n", main="Protein intensity per number of constructs bound")

##
## estudi d'homologia en humans
##
# baixar els homolegs en human dels binders
hom <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene"),
             filters="ensembl_gene_id", values=unique(unlist(strsplit(ids.m, ";"))), 
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))
ann.hs <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "name_1006"), # name_1006 is the GO_term_description
                filters="ensembl_gene_id", values=hom$hsapiens_homolog_ensembl_gene,
                mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="hsapiens_gene_ensembl"))
hom$hsapiens_homolog_gene_name <- sapply(hom$hsapiens_homolog_ensembl_gene, function(x) {
    ann.hs$external_gene_name[match(x, ann.hs$ensembl_gene_id)]
})
hom$hsapiens_homolog_goslim    <- sapply(hom$hsapiens_homolog_ensembl_gene, function(x) {
    go <- ann.hs$name_1006[match(x, ann.hs$ensembl_gene_id)]
    paste(go, collapse=";")
})

write.csv(hom, file=paste0(OUT, "/heatmap_1+1_threshold.hom.csv"))

# fer un parell de plots amb els GO terms
p <- ggplot(as.data.frame(table(hom$hsapiens_homolog_goslim[hom$hsapiens_homolog_goslim != "NA"])), aes(x=Var1, y=Freq)) +
        geom_bar(stat="identity") + coord_flip() + xlab("") + theme_bw()
p

library(wordcloud)
wordcloud(hom$hsapiens_homolog_goslim, colors=brewer.pal(6,"Dark2"), main="human homologs GO terms")

dev.off()
