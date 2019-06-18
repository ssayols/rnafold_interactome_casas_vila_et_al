############################
##
## Read FACS data and do a couple of plots:
##   #-heatmap table, all versus all GFP ratios
##   -ordered scatter of log2 enrichment ratio vs. GFP
##   -summary boxplot of log2 enrichment ratio vs. GFP, all constructs
##   -individual boxplot of log2 enrichment ratio vs. GFP
##
############################
library(readxl)

#pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal <- paste0(RColorBrewer::brewer.pal(8, "Dark2")[2:1], "AA")
#BASEDIR <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/facs"
#setwd(BASEDIR)
BASEDIR <- getwd()

# get command args
#args <- commandArgs(trailingOnly=TRUE)
args <- c("./20171017/17-1017FlowJo_plate4_Znshock.xls", "./20171017/17-1017FlowJo_plate4_Znshock")
if(length(args) != 2 || !file.exists(args[1]))
    stop("Call with: Rscript facs.R ./facs/FlowJoAnalysis.xls ./facs/FlowJoAnalysis_out_prefix\n")

pdf(paste0(args[2], ".pdf"))

##
## read a sheet
##
x <- as.data.frame(read_xls(args[1], 1))
fdata <- x$`yeast/Singlets/living | Median (BL488nm 530_30-A)`
names(fdata) <- gsub("(^Specimen_001_|_[0-9]+\\.fcs$)", "", x[, 1])
names(fdata) <- gsub(" ", "_", names(fdata))
names(fdata) <- gsub("NULL", "Null", names(fdata))

# enrichment over GFP; check individual replicate enrichment over avg. GFP
gfp <- grepl("^GFP_", names(fdata))
samples <- sub("^([A-Za-z0-9]+)_.+", "\\1", names(fdata))
ratios    <- sapply(fdata, function(x) x / mean(fdata[gfp]))
ratio.avg <- tapply(fdata, samples, function(x) mean(x / mean(fdata[gfp])))
ratio.sd  <- tapply(fdata, samples, function(x) sd(x / mean(fdata[gfp])))
ratio.sd  <- ifelse(is.na(ratio.sd), 0, ratio.sd)
ratio.se  <- tapply(fdata, samples, function(x) sd(x / mean(fdata[gfp])) / sqrt(length(x)))
ratio.se  <- ifelse(is.na(ratio.se), 0, ratio.se)
ratio.val <- tapply(fdata, samples, function(x) paste(x / mean(fdata[gfp]), collapse=":"))

##
##   -heatmap table, all versus all GFP ratios
##
#x <- tapply(fdata, samples, mean)
#ratio.all.vs.all <- sapply(x, function(y) log2(x / y))
#df <- reshape2::melt(ratio.all.vs.all)
#
#ggplot(df, aes(x=Var1, y=Var2, size=value)) +
#  geom_point(aes(colour=value)) +
#  scale_colour_distiller(direction=-1) +
#  xlab("") + ylab("") + theme_bw() +
#  theme(axis.text.x=element_text(angle=90, hjust=1))

##
##   -ordered scatter of log2 enrichment ratio vs. GFP
##
# reorder
o <- order(ratio.avg)
ratio.avg <- ratio.avg[o]
ratio.sd  <- ratio.sd [o]
ratio.se  <- ratio.se [o]
ratio.val <- ratio.val[o]

# and plot
plot(1:length(ratio.avg), log2(ratio.avg), pch=16, cex=.6, col="#00000050", xaxt="n",
     xlab="", ylab="log2 ratio construct / GFP")
abline(h=0, lty=2, col="blue")
text(1:length(ratio.avg), log2(ratio.avg), labels=names(ratio.avg), pos=3, srt=90, cex=.6)

##
##   -summary boxplot of log2 enrichment ratio vs. GFP, all constructs
##
do.barplot <- function(x, se, ...) {
    b <- barplot(x, ylim=range(c(0, x, x - se * 2, x + se * 2)), ...)
    segments(b, x - se * 2, b, x + se * 2, lwd=1.5)
    suppressWarnings(arrows(b, x - se * 2, b, x + se * 2, angle=90, code=3, length=0.05, lwd=1.5))
    box(bty="l")
}

# and plot
do.barplot(ratio.avg, ratio.se, las=2, cex.axis=.6, cex.names=.6,
           col=pal[as.numeric(!grepl("^GFP$", names(ratio.avg))) + 1],
           ylab="ratio construct / GFP")
do.barplot(log2(ratio.avg), ratio.se, las=2, cex=.6, cex.names=.6,
           col=pal[as.numeric(!grepl("^GFP$", names(ratio.avg))) + 1],
           ylab="log2 ratio construct / GFP")

##
##   -individual boxplot of log2 enrichment ratio vs. GFP
##
par(mfrow=c(2, 2))
invisible(
    Map(function(x, se, construct) {
        gfp <- grepl("^GFP$", names(ratio.avg))
        x <- c(ratio.avg[gfp], x)
        names(x) <- c("GFP", construct)
        se <- c(ratio.se[gfp], se)

        do.barplot(x, se, las=2, space=.5, col=pal[1:2], ylab="ratio construct / GFP")

        replicates <- grep(paste0("^", construct, "_"), names(ratios))
        if(length(replicates) > 2) {
          legend("topright", bty="n", cex=3/4,
                 legend=paste("pval=",
                              format.pval(t.test(ratios[grep("^GFP_", names(ratios))],
                                                 ratios[replicates])$p.value)))
        }
    }, ratio.avg, ratio.se, names(ratio.avg))
)

write.csv(data.frame(ratio.avg=ratio.avg,
                     ratio.se =ratio.se,
                     ratio.sd =ratio.sd,
                     ratio.val=ratio.val),
          file=paste0(args[2], ".csv"))

dev.off()
