################################
##
## plot number of binders per construct type (3', CDS, 5')
## this goes to one of the pannels of figure 1
##
###############################
library(ggplot2)
library(reshape)

WD <- "/fsimb/groups/imb-buttergr/=EVOREG/rnalib"
setwd(WD)

##
## read heatmap matrix with constructs and binders
##
# read heatmap of binders and constructs
h <- read.csv("./analysis/heatmap_1+1_threshold.csv")
h <- h[-1:-4, ] # remove controls
h <- do.call(rbind, apply(h, 1, function(x) do.call(rbind, lapply(0:sum(grepl(";", x[1])), function(i) x)))) # duplicate protein groups
h <- as.matrix(h[, -1:-2])
mode(h) <- "double"
#h <- h[, apply(h, 2, function(x) any(x > 0))] # remove constructs with no binders

##
## plot number of binders per construct type (3', CDS, 5')
## 
# sum binders per construct and add position
df <- melt(apply(h, 2, function(x) sum(x > 0)))
consts.dict  <- read.csv("./constructs/constructs.csv")
df$pos <- consts.dict$pos[match(rownames(df), consts.dict$name)]  # some are duplicated

# add rank within pos
df <- do.call(rbind, lapply(split(df, df$pos), function(x) {
  x$rank <- rank(x$value, ties="first")
  x$rank_norm <- x$rank / max(x$rank)   # normalize rank between 0..1
  x
}))

pdf("./analysis/plot_binders.pdf", width=10, height=10)

lapply(list(all=df, with_binders=subset(df, value > 0)), function(df) {
  print(
    ggplot(df, aes(x=rank_norm, y=value, color=pos)) +
      geom_point(size=3, alpha=.8) +
      scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
      theme_bw()
  )

#  print(
#    ggplot(df, aes(x=value, fill=pos)) +
#      geom_density(alpha=.8) +
#      scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
#      theme_bw()
#  )

  print(
    ggplot(df, aes(x=rank, y=value)) +
      geom_point(size=3, alpha=.8) +
      theme_bw() +
      facet_wrap(~pos, scales="free_x", ncol=3)
  )

  print(
    ggplot(df, aes(x=rank, y=value)) +
      geom_point(size=3, alpha=.8) +
      theme_bw() +
      facet_wrap(~pos, scales="free_x", nrow=3)
  )

  invisible(0)
})

dev.off()
