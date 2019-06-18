################################
##
## Integration of:
##   -the CellMap data (http://thecellmap.org/costanzo2016/) and
##   #-String PPI DB (http://string-db.org/cgi/download.pl)
##   #-BioGrid (integrates both string and the CellMap
##
## Take our RBP-->construct interaction data and:
##   -add the genetic interaction info from the CellMap on our RBP-->construct interaction table
##   -GO enriched terms for all CellMap interactors of our RBP
##   -GO enriched terms of similar RBP profiles in the CellMap
##
################################
library(parallel)
#library(plotly)
library(gplots)
library(ggplot2)
library(reshape)
library(WriteXLS)
library(igraph)
library(RColorBrewer)

CORES <- 16
SCORE <- .08 # as described as |E| (genetic interaction score) in http://science.sciencemag.org/content/327/5964/425.full
PCC   <- .2  # as described before: the Pearson Correlation Coefficient between genes with similar genetic interaction profiles
WD <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib"
setwd(WD)

pdf("./analysis/cellmap.pdf", width=10, height=10)
pal  <- colorRampPalette(c("white", "black"))(100)
pal2 <- c("#E69F00", "#56B4E9", "#009E73")

##
## read heatmap matrix with constructs and binders
##
# read dictionaries
consts.dict  <- read.csv("./constructs/constructs.csv")
ensembl.dict <- read.delim("./functional_follow_up/ensembl_dic.txt")
goslim       <- read.delim("./analysis/db/go_slim_mapping.tab", comment.char="#")
goslim       <- goslim[goslim$GOtype == "P", ]  # only BP
goslim       <- goslim[goslim$term != "biological_process", ]

# read heatmap of binders and constructs
h <- read.csv("./analysis/heatmap_1+1_threshold.csv")
h <- h[-1:-4, ] # remove controls
binders <- h[, 1]
h <- do.call(rbind, apply(h, 1, function(x) do.call(rbind, lapply(0:sum(grepl(";", x[1])), function(i) x)))) # duplicate protein groups
binders <- unlist(strsplit(binders, ";"))
h <- as.matrix(h[, -1:-2])
mode(h) <- "double"
h <- h[, apply(h, 2, function(x) any(x > 0))] # remove constructs with no binders
constructs <- consts.dict$gene[match(colnames(h), consts.dict$name)]

# read external data
cellmap <- as.matrix(read.csv("./analysis/db/thecellmap.csv", row.names=1))
#string  <- as.matrix(read.csv("./analysis/db/string.csv", row.names=1))
#biogrid <- as.matrix(read.csv("./analysis/db/biogrid.csv", row.names=1))

##
##   -add the genetic interaction info from the CellMap on our RBP-->construct interaction table
##
# select from cellmap the bindersXconstructs
cellmap.h <- cellmap[rownames(cellmap) %in% binders, colnames(cellmap) %in% constructs]
h.cellmap <- h[match(rownames(cellmap.h), binders), match(colnames(cellmap.h), constructs)]
h.cellmap <- ifelse(h.cellmap > 0, 1, 0)    # convert our interaction matrix into qualitative values
rownames(h.cellmap) <- rownames(cellmap.h)
colnames(h.cellmap) <- colnames(cellmap.h)

# remove binders which don't have any information left about their bound constructs in the Cellmap
out <- apply(h.cellmap, 1, function(x) all(x == 0))
cellmap.h <- cellmap.h[!out, ]
h.cellmap <- h.cellmap[!out, ]

# calculate row/col order in the heatmap based on hierarchical clustering
cluster.rows <- hclust(dist(h.cellmap))
cluster.cols <- hclust(dist(t(h.cellmap)))

#plot_ly() %>% add_surface(z           =h.cellmap[cluster.rows$order, cluster.cols$order],
#                          surfacecolor=abs(cellmap.h[cluster.rows$order, cluster.cols$order]))

heatmap.2(abs(cellmap.h[cluster.rows$order, cluster.cols$order]), trace="none", scale="none", col=pal, Colv=NA, Rowv=NA, dendrogram="none", main="thecellmap")
heatmap.2(h.cellmap[cluster.rows$order, cluster.cols$order], trace="none", scale="none", col=pal, Colv=NA, Rowv=NA, dendrogram="none", main="binders")

# identify on the h.cellmap which cells are also enriched in the cellmap.h dataset
h.cellmap.melt <- melt.matrix(h.cellmap)
cellmap.h.melt <- melt.matrix(cellmap.h)

annotate <- function(x) paste(x, ensembl.dict$Associated.Gene.Name[match(x, ensembl.dict$Ensembl.Gene.ID)], sep=" - ")

df <- data.frame(x=factor(annotate(cellmap.h.melt$X1), levels=annotate(rownames(h.cellmap)[cluster.rows$order])),
                 y=factor(annotate(cellmap.h.melt$X2), levels=annotate(colnames(h.cellmap)[cluster.cols$order])),
                 value=cellmap.h.melt$value * h.cellmap.melt$value)  # our interaction table is TRUE/FALSE (1/0)
df <- df[!is.na(df$value), ]
df <- df[abs(df$value) > 0, ]
df$score <- factor(ifelse(df$value >  SCORE, "positive",
                   ifelse(df$value < -SCORE, "negative", "below threshold")),
                   levels=c("positive", "negative", "below threshold"))

p <- ggplot(df, aes(x=x, y=y, color=score, size=abs(value))) +
       geom_point() +
       scale_x_discrete("binders") +
       scale_y_discrete("constructs") +
       scale_color_manual("interaction score", values=c("green", "red", "grey")) +
       scale_size_continuous("", range=c(0,10), limits=c(0, 1)) +
       theme_bw() +
       theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
print(p)
print(p + guides(color=FALSE, size=FALSE) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()))

##
##   -GO enriched terms for all CellMap interactors of our RBP
##
# make cluster
cl <- makeCluster(CORES)
x <- clusterEvalQ(cl, library(org.Sc.sgd.db))
x <- clusterEvalQ(cl, library(GOstats))
clusterExport(cl, c("binders", "cellmap", "ensembl.dict", "SCORE"))

# calculate the GO enriched terms of the genetic interactors for our binders
GO.genetic.interactors <- parLapply(cl, binders[binders %in% rownames(cellmap)], function(binder) {

    # get from our lists of binders which have
    genetic.interactors <- colnames(cellmap)[!is.na(cellmap[binder, ]) & abs(cellmap[binder, ]) > SCORE]

    # calcualte hyperG test
    param2  <- new("GOHyperGParams", geneIds=genetic.interactors, universeGeneIds=Lkeys(org.Sc.sgdGO), annotation="org.Sc.sgd.db", ontology="BP")
    hyp     <- hyperGTest(param2)
    result  <- summary(hyp, categorySize=2)
    result  <- data.frame(result, FDR=p.adjust(result$Pvalue, "fdr"))
    result$genes <- sapply(geneIdUniverse(hyp)[result$GOBPID], function(x) paste(unique(genetic.interactors[genetic.interactors %in% x]), collapse="; "))
    result$genes_names <- sapply(geneIdUniverse(hyp)[result$GOBPID], function(x)
        paste(unique(ensembl.dict$Associated.Gene.Name[match(genetic.interactors[genetic.interactors %in% x], ensembl.dict$Ensembl.Gene.ID, nomatch=0)]), collapse="; "))

    result
})

names(GO.genetic.interactors) <- binders[binders %in% rownames(cellmap)]
GO.genetic.interactors <- GO.genetic.interactors[order(names(GO.genetic.interactors))]
WriteXLS("GO.genetic.interactors", "./analysis/cellmap.GO_interactors.xlsx")

##
##   -GO enriched terms of similar RBP profiles in the CellMap
##
ensembl.dict$Associated.Gene.Name <- ifelse(ensembl.dict$Associated.Gene.Name == "", ensembl.dict$Ensembl.Gene.ID, ensembl.dict$Associated.Gene.Name)

mapply(function(genes, f, color) {
  set.seed(666)

  # calculate similarities between interaction profiles in our RBP
  pcc <- cor(t(cellmap[rownames(cellmap) %in% genes, ]), use="pairwise.complete.obs")
  pcc <- ifelse(is.na(pcc) | abs(pcc) > .8, 0, pcc)    # remove artifacts due to not enough data
  for(i in 1:nrow(pcc)) pcc[i, i] <- 1  # bring the diagonal back to 1

  # heatmap with the correlation, in the meantime: pdf("./analysis/cellmap_dist_heatmap.pdf", width=10, height=10)
  x <- 1 - abs(pcc)
  ccol <- ifelse(colnames(x) %in% binders & colnames(x) %in% constructs, pal2[3],
          ifelse(colnames(x) %in% binders, pal2[1], pal2[2]))
  heatmap.2(x, col=pal, trace="none", ColSideColors=ccol, RowSideColors=ccol,
            distfun=function(x) dist(x, method="euclidean"),
            hclustfun=function(d) hclust(d, method="complete"))

  # build network and detect communities
  pcc.binary <- ifelse(abs(pcc) > PCC, 1, 0)  # only link nodes with an score above the threshold
  for(i in 1:nrow(pcc.binary)) pcc.binary[i, i] <- 0        # remove autoreferences
  out <- apply(pcc.binary, 1, function(x) all(x == 0))
  pcc.binary <- pcc.binary[!out, !out]
  colnames(pcc.binary) <- rownames(pcc.binary) <- annotate(colnames(pcc.binary))
  net <- graph_from_adjacency_matrix(pcc.binary, mode="undirected")
  clp <- cluster_label_prop(net)
  V(net)$community <- clp$membership
  V(net)$color     <- if(color) {
                        nodes <- gsub(" - .*$", "", names(V(net)))
                        ifelse(nodes %in% binders & nodes %in% constructs, pal2[3],
                        ifelse(nodes %in% binders, pal2[1], pal2[2]))
                      } else {
                        pal2[1]
                      }
  # color edges in red if there's interaction between the protein and the rna; black otherwise
  E(net)$color     <- if(color) {
                        nodes <- sapply(E(net), function(x) V(net)[inc(x)]$name)  # get the A--B pair of nodes for each edge
                        node1 <- gsub(" - .*$", "", nodes[1, ])
                        node2 <- gsub(" - .*$", "", nodes[2, ])
                        mapply(node1, node2,FUN=function(n1, n2) {
                          if(n1 %in% binders && n2 %in% constructs) {
                            if(any(h[binders == n1, constructs == n2] > 0)) "red" else "black"  # there's some duplicated constructs
                          } else if(n2 %in% binders && n1 %in% constructs) {
                            if(any(h[binders == n2, constructs == n1] > 0)) "red" else "black"  # there's some duplicated constructs
                          } else {
                            "black"
                          }
                        })
                      } else {  # don't color
                        "black"
                      }

  invisible(
    lapply(list(set.vertex.attribute(net, "label", value=gsub(" - .*$", "", names(V(net)))),   # display the gene ID
                set.vertex.attribute(net, "label", value=gsub("^.+ - ", "", names(V(net))))),  # display the gene name
           function(net) {
      set.seed(666)
      plot(clp, net, vertex.frame.color="white", col=V(net)$color, mark.col="#00000010", mark.border="#000000AA",
           vertex.size=10, vertex.label.family="Helvetica", vertex.label.cex=2/4, vertex.label.color="black",
           edge.color=E(net)$color,
           layout=layout_nicely)#with_kk) # Kamada-Kawai network layout (aka the spring-embedded network layout)

      if(color) legend("bottomright", c("RBP only", "mRNA only", "RBP+mRNA"), fill=pal2, border=NULL, cex=3/4, bty="n")
    })
  )

  # calculate GO on the communities detected
  x <- clusterEvalQ(cl, rm(clp))
  clusterExport(cl, "clp", envir=environment())
  GO.communities <- parLapply(cl, unique(membership(clp)), function(i) {

      # get from our lists of genes which have
      genes <- strsplit(clp$names[clp$membership == i], " - ")
      genes <- sapply(genes, function(x) x[1])

      # calcualte hyperG test
      param2  <- new("GOHyperGParams", geneIds=genes, universeGeneIds=Lkeys(org.Sc.sgdGO), annotation="org.Sc.sgd.db", ontology="BP")
      hyp     <- hyperGTest(param2)
      result  <- summary(hyp, categorySize=2)
      result  <- data.frame(result, FDR=p.adjust(result$Pvalue, "fdr"))
      result$genes <- sapply(geneIdUniverse(hyp)[result$GOBPID], function(x) paste(unique(genes[genes %in% x]), collapse="; "))
      result$genes_names <- sapply(geneIdUniverse(hyp)[result$GOBPID], function(x)
          paste(unique(ensembl.dict$Associated.Gene.Name[match(genes[genes %in% x], ensembl.dict$Ensembl.Gene.ID, nomatch=0)]), collapse="; "))

      result
  })

  names(GO.communities) <- paste("community", unique(membership(clp)))
  GO.communities <- GO.communities[order(names(GO.communities))]
  WriteXLS("GO.communities", f)
}, list(binders, constructs, c(binders, constructs)),
   list("./analysis/cellmap.GO_RBP_communities.xlsx", "./analysis/cellmap.GO_constructs_communities.xlsx", "./analysis/cellmap.GO_RBP+constructs_communities.xlsx"),
   color=c(FALSE, FALSE, TRUE))

stopCluster(cl)
dev.off()
