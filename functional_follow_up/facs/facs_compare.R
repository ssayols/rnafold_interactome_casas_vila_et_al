###############################
##
## Compare enrichment of constructs vs. GFP between plates (growing conditions)
## As input, take the results from facs.R
##
###############################
pal <- paste0(RColorBrewer::brewer.pal(8, "Dark2")[c(2, 1, 3:8)], "AA")
#BASEDIR <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up"
#setwd(BASEDIR)
BASEDIR <- getwd()

# get command args
runstr <- "Call with: Rscript facs_compare.R ./facs/ctl.csv ./facs/plate2.csv [./facs/plate3.csv ...] ./facs/zinc_out_suffix\n"
#args <- commandArgs(trailingOnly=TRUE)
args <- c("./20171017/20171013FlowJo_plate1.csv",
          "./20171017/17-1017FlowJo_Plate1_SD-G.csv",
          "./20171017/17-1017FlowJo_Plate2_SD-N.csv",
          "./20171017/17-1017FlowJo_plate3_LZM.csv",
          "./20171017/17-1017FlowJo_plate4_Znshock.csv",
          "./20171017/17-1017_cross_comparison")

if(length(args) < 3)
  stop(runstr)

OUT <- args[length(args)]
EXPERIMENTS <- args[-length(args)]
files.exist <- sapply(EXPERIMENTS, file.exists)
if(!all(files.exist))
  stop("some of the files don't exist: ", args[!files.exist])

pdf(paste0(OUT, ".pdf"))

##
## read sheets
##
plates <- gsub("(^[0-9\\-]+FlowJo_|\\.csv$)", "", basename(EXPERIMENTS))
ratios <- Reduce(function(x, y) merge(x, y, by=1, all=F), lapply(EXPERIMENTS, read.csv))
colnames(ratios) <- c("construct", paste(c("ratio.avg", "ratio.se", "ratio.sd", "ratio.val"),
                                         rep(plates, each=4)))

##
## barplot ratios
##
do.barplot <- function(x, se, ...) {
  b <- barplot(x, ylim=range(c(0, x, x - se * 2, x + se * 2)), ...)
  segments(b, x - se * 2, b, x + se * 2, lwd=1.5)
  suppressWarnings(arrows(b, x - se * 2, b, x + se * 2, angle=90, code=3, length=0.05, lwd=1.5))
  box(bty="l")
}

# the plot
par(mfrow=c(2, 2))
invisible(
  apply(ratios, 1, function(x) {
    # do the plot
    avg <- c(1, as.numeric(x[grep("^ratio.avg", names(x))]))
    names(avg) <- c("GFP", plates)
    do.barplot(x =avg,
               se=c(0, as.numeric(x[grep("^ratio.se" , names(x))])),
               main=x["construct"],
               las=2, space=.5, col=pal[1:(length(plates) + 1)],
               ylab="ratio construct / GFP")

    # add p-value
    y <- reshape2::melt(lapply(x[grep("^ratio.val", names(x))], function(x) as.numeric(unlist(strsplit(x, ":")))))
    legend("topright", bty="n", cex=3/4,
           legend=paste("pval=", format.pval(anova(lm(y$value ~ y$L1))[1, "Pr(>F)"])))
  })
)

write.csv(ratios, row.names=FALSE,
          file=paste0(OUT, ".csv"))

dev.off()

