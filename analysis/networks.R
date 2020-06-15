################################
##
## infer communities and GO BP associated to such communities
## from the network analysis of the constructs and their binders
##
################################
library(igraph)
library(RColorBrewer)

WD <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib"
setwd(WD)

pdf("./analysis/network.pdf")

##
## build the vertex/edges matrices
##
consts.dict  <- read.csv("./constructs/constructs.csv")
ensembl.dict <- read.delim("./functional_follow_up/ensembl_dic.txt", sep="\t")
h <- read.csv("./analysis/heatmap_1+1_threshold.csv")
h <- h[-1:-4, ] # remove controls
rownames(h) <- h$X
binders     <- unlist(strsplit(h$X, ";"))
h <- as.matrix(h[, -1:-2])
h <- h[, apply(h, 2, function(x) any(x > 1))]
constructs  <- consts.dict$gene[match(colnames(h), consts.dict$name)]

# vertex
vertex <- data.frame(id=unique(c(binders, constructs)))
vertex$type <- ifelse(vertex$id %in% constructs & vertex$id %in% binders, "both",
               ifelse(vertex$id %in% constructs, "construct", "binder"))
vertex <- do.call(rbind, by(vertex, vertex$id, function(x) data.frame(id=unlist(strsplit(x$id, ";")), type=x$type)))
vertex$gene <- ensembl.dict$Associated.Gene.Name[match(vertex$id, ensembl.dict$Ensembl.Gene.ID)]
vertex$gene <- ifelse(is.na(vertex$gene) | vertex$gene == "", vertex$id, vertex$gene)

# edges
edges.wells <- apply(h, 1, function(x) colnames(h)[x > 1])
edges.genes <- apply(h, 1, function(x) constructs [x > 1])
out   <- sapply(edges.wells, length) == 0  # remove controls
edges <- mapply(function(from, to.genes, to.wells) { data.frame(from=from, to=to.genes, weight=h[from, to.wells]) },
                rownames(h)[!out], edges.genes[!out], edges.wells[!out], SIMPLIFY=FALSE)
edges <- do.call(rbind, edges)
edges <- do.call(rbind, by(edges, edges$from, function(x) data.frame(from=unlist(strsplit(x$from, ";")), to=x$to, weight=x$weight)))

##
## draw network
##
net <- graph_from_data_frame(d=edges, vertices=vertex, directed=T)
E(net)$width <- 1 + (E(net)$weight - min(E(net)$weight)) / (max(E(net)$weight) - min(E(net)$weight))    # scale weights between [0,1]
l <- layout_with_fr(net)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

# Community detection based on label propagation:
clp <- cluster_walktrap(net, weights=E(net)$weight)
V(net)$community <- clp$membership

# plot
pal <- paste0(brewer.pal(9, "Set1"), "80")
plot(net, xlim=c(-2, 2), ylim=c(-2, 2), edge.arrow.size=.2, edge.color="orange", rescale=F,
     vertex.color=pal[factor(V(net)$community)], vertex.frame.color="white", vertex.label="",
     edge.curved=.3, layout=l*2)
#     vertex.label=V(net)$gene, vertex.label.color=pal[as.factor(V(net)$type)])

legend("bottomleft", legend=levels(factor(V(net)$type)), fill=pal[factor(levels(factor(V(net)$type)))], pch=21, bty="n")

dev.off()
