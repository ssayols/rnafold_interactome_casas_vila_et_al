################################
##
## Calculate the in-vivo vs. in-vitro similarity of the structures
##
## It's the relationship between the amount of structure, computed as
## the Gini coefficient, and the correlation between the in-vivo vs. in-vitro
## computed as the R-squared coefficient of the Spearman correlation.
##
## The raw data is normalized as described in the article: proportionally
## to the most reactive base within the given structure. With a small
## modification: no 50-200bp windows, we take the complete structure.
##
################################
options(stringsAsFactors=F)
set.seed(666)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
#library(ggrepel)
library(RColorBrewer)
library(reshape)

NBOOT <- 200
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
pdf("./DMS-seq_data/vivo_vs_vitro.pdf")

# compute gini coefficient as in: http://mathworld.wolfram.com/GiniCoefficient.html
gini <- function(x) {
    sum(sapply(1:length(x), function(i) sum(abs(x[i] - x)))) / (2 * length(x)^2 * mean(x))
}

##
## read input data and create GRanges objects
##
# read and parse DMS data. Format is wiggle
dms.vivo.minus  <- import.wig("./DMS-seq_data/Silvi/May2013_VivoAllextra_minus.txt")
dms.vivo.plus   <- import.wig("./DMS-seq_data/Silvi/May2013_VivoAllextra_plus.txt")
dms.vitro.minus <- import.wig("./DMS-seq_data/Silvi/May2013_VitroAllextra_minus.txt")
dms.vitro.plus  <- import.wig("./DMS-seq_data/Silvi/May2013_VitroAllextra_plus.txt")
dms.denat.minus <- import.wig("./DMS-seq_data/Silvi/May2013_DenaturedAllextra_minus.txt")
dms.denat.plus  <- import.wig("./DMS-seq_data/Silvi/May2013_DenaturedAllextra_plus.txt")

# read input structures
structs <- read.delim("./primers/feat.csv", head=F)
structs <- with(structs, GRanges(seqnames=V1, ranges=IRanges(V2, V3), strand=V5, name=V4))

# get NBOOT random regions with dms signal
structs.random <- replicate(NBOOT, {
    strands <- table(strand(structs))
    s <- rbinom(1, 1, strands[1] / sum(strands))
    x <- sample(switch(s+1, dms.vivo.minus, dms.vivo.plus), 1)
    strand(x) <- if(s == 1) "+" else "-"
    end(x) <- end(x) + rpois(1, mean(width(structs)))
    x
})
structs.random <- unlist(GRangesList(structs.random))

##
## Compute the in vivo/vitro data for each structure
##
get.dms.data <- function(x) {
    # get in vivo/vitro data
    dms.vivo  <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.vivo.minus, x),
                                      "+"=subsetByOverlaps(dms.vivo.plus , x)))
    dms.vitro <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.vitro.minus, x),
                                      "+"=subsetByOverlaps(dms.vitro.plus , x)))
    dms.denat <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.denat.minus, x),
                                      "+"=subsetByOverlaps(dms.denat.plus , x)))

    # normalize to the most reactive base
    dms.vivo$score  <- dms.vivo$score  / max(dms.vivo$score)
    dms.vitro$score <- dms.vitro$score / max(dms.vitro$score)
    dms.denat$score <- dms.denat$score / max(dms.denat$score)

    # merge by pos (dropping incomplete positions) and calculate r-squared and gini coefficient
    dms <- Reduce(function(x, y) merge(x, y, by=c("start", "end")), list(dms.vivo, dms.vitro, dms.denat))
    r   <- cor(dms$score.x - dms$score, dms$score.y - dms$score, method="spearman")
    g.vivo  <- gini(dms$score.x)
    g.vitro <- gini(dms$score.y)
    g.denat <- gini(dms$score)

    c(r=r, g.vivo=g.vivo, g.vitro=g.vitro, g.denat=g.denat)
}

dms <- as.data.frame(t(sapply(structs, get.dms.data)))
dms.random <- as.data.frame(t(sapply(structs.random, get.dms.data)))
rownames(dms) <- structs$name

##
## plot the results
##
# in vivo/vitro correlation
p <- ggplot(dms, aes(x=g.vivo - g.vitro, y=r^2)) +
        geom_point(size=3, alpha=1/3) +
        xlab("Gini difference (in vivo - in vitro)") +
        ylab("R-squared") +
        ggtitle("Correlation of in vivo/vitro structures") +
        theme_bw()
print(p)

# in vivo/vitro correlation, structure dependant
p <- ggplot(dms, aes(x=g.vivo - g.denat, y=g.vitro - g.denat)) +
        geom_point(aes(colour=r^2), size=3, alpha=2/3) +
        geom_smooth(method="lm") +
        scale_colour_distiller("R-squared", palette="Spectral") +
        xlab("Gini difference (in vivo - denatured)") +
        ylab("Gini difference (in vitro - denatured)") +
        ggtitle("Correlation of in vivo/vitro structures") +
        theme_bw()
print(p)

# same plot, but adding random in a second facet
dms$type <- "Weissman"
dms.random$type <- "random"
x <- rbind(dms, dms.random)
x$r <- x$r^2
p <- ggplot(x, aes(x=g.vivo - g.denat, y=g.vitro - g.denat)) +
        geom_point(aes(colour=r^2), size=3, alpha=2/3) +
#        geom_smooth(method="lm") +
        scale_colour_distiller("R-squared", palette="RdBu", direction=-1) +
        xlab("Gini difference (in vivo - denatured)") +
        ylab("Gini difference (in vitro - denatured)") +
        ggtitle("Correlation of in vivo/vitro structures") +
        theme_bw() +
        facet_wrap(~type)
print(p)

# in vivo/vitro correlation, structure dependant with random
dms$type <- "Weissman"
dms.random$type <- "random"
x <- rbind(dms, dms.random)
x$r <- x$r^2

p <- ggplot(x, aes(x=g.vivo - g.denat, y=g.vitro - g.denat, colour=type)) +
        geom_point(aes(size=r)) +
#        geom_smooth(method="lm") +
        scale_colour_manual("", values=c("#56B4E950", "#E69F0050")) +
        scale_size_continuous("R-squared", range=c(0, 8)) +
        xlab("Gini difference (in vivo - denatured)") +
        ylab("Gini difference (in vitro - denatured)") +
        ggtitle("Correlation of in vivo/vitro structures") +
        theme_bw()
print(p)

# same plot, using base R graphics
plot(x=x$g.vivo  - x$g.denat,
     y=x$g.vitro - x$g.denat,
     xlab="Gini difference (in vivo - denatured)",
     ylab="Gini difference (in vitro - denatured)",
     pch=16,
     col=ifelse(x$type == "Weissman", alpha("orangered", .8), alpha("black", .8)))

# show that weissman's structures have higher in-vivo - in-vitro correlation
#boxplot(outline=FALSE, lty=1, staplewex=0, boxwex=0.8, boxlwd=1, medlwd=1, col=cols, add=TRUE, xaxt="n", yaxt="n")
boxplot(x$r ~ x$type,
        main="in-vivo ~ in-vitro correlation",
        ylab="R^2",
        ylim=c(0, 1),
        col=alpha("grey", .8), lty=1, boxlwd=4)

# validation of the data: capacity of denatured/in vivo/in vitro to form structures
x <- melt.data.frame(x[, -1]) # remove Rho and keep only gini's
x$variable <- ifelse(x$variable == "g.vivo",  "in vivo", 
              ifelse(x$variable == "g.vitro", "in vitro", "denatured"))

f <- function(x) {
    p <- ggplot(x, aes(x=value)) +
            geom_density(aes(fill=variable), alpha=1/3) +
            scale_fill_manual("Experiment", values=c("#E69F00", "#56B4E9", "#009E73")) +
            ggtitle("Capacity to form structures") +
            xlab("Gini coefficient") +
            theme_bw()
}
print(f(x[x$type == "Weissman",]))
print(f(x[x$type == "random",]))

dev.off()
