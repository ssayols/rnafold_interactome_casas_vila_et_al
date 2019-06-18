##
## plot the result of the GO term slim analysis of the constructs
##
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape)
library(RColorBrewer)

TOP=25

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/GO/constructs")
pdf("GOterms.constructs.pdf")

# primer fer un barplot amb la distribucio de constructes per 5'UTR, CDS, 3'UTR
constructs <- read.csv("../../../constructs/constructs.csv")
constructs$gene <- gsub("_\\d+_(CDS|UTR5|UTR3)_(S|L)$", "", constructs$name1)
constructs$pos  <- gsub("^.+(CDS|UTR5|UTR3).+", "\\1", constructs$name1)
y <- table(constructs$pos)[c("UTR5", "CDS", "UTR3")]
x <- barplot(y, space=.6, width=.4, col=RColorBrewer::brewer.pal(3, "Set2"), cex.names=1.5, cex.axis=1.5, ylim=c(0, max(y) + 20))
text(x, y + 10, y, cex=1.5)

# llegir el resultat del GO slim mapper de SGD amb les proteines dels constructes
goslim <- read.delim("slimMapper_SGD_constructs.txt")
goslim$observed      <- as.numeric(gsub("(\\d+) out of (\\d+) genes,.+","\\1",goslim$Frequency))
goslim$totalObserved <- as.numeric(gsub("(\\d+) out of (\\d+) genes,.+","\\2",goslim$Frequency))
goslim$expected      <- as.numeric(gsub("(\\d+) of (\\d+) genes,.+","\\1",goslim$Genome.Frequency))
goslim$totalExpected <- as.numeric(gsub("(\\d+) of (\\d+) genes,.+","\\2",goslim$Genome.Frequency))
goslim$ratio         <- (goslim$observed / goslim$totalObserved) / (goslim$expected / goslim$totalExpected)

filt   <- grepl("other", goslim$GO.term) | goslim$expected == 0
goslim <- goslim[!filt, ]

p <- ggplot(goslim[rev(order(goslim$ratio))[1:TOP], ], aes(x=GO.term, y=ratio)) + 
     geom_bar(stat="identity", position="dodge") +
     scale_fill_brewer(palette="Set2") +
     scale_y_continuous("log2 O/E ratio") +
     scale_x_discrete("") +
     coord_flip() + theme_bw() +
     theme(text=element_text(size=7), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7))
p

dev.off()
