library(ggplot2)
library(gridExtra)

x <- read.delim("raw.txt", comment.char="#")
# R> head(data)
#   binder loop condition type sample input       Cq
# 1   Pab1 RBG2  wildtype  +RT   NC13     0 30.88703
# 2   Pab1 RBG2  wildtype  +RT   NC13     0 30.23182
# 3   Pab1 RBG2  wildtype  +RT   NC13     0 29.96186
# 4   Pab1 YEF3  wildtype  +RT   NC13     0 26.95399
# 5   Pab1 YEF3  wildtype  +RT   NC13     0 25.95151
# 6   Pab1 YEF3  wildtype  +RT   NC13     0 25.60303

#pdf("rip.pdf")
#par(mfrow=c(1, 4))

p <- lapply(split(x, x$loop), function(x) {
    lapply(split(x, x$binder), function(x) {
        # calculate percentage input of +RT & wildtype condition
        input_rt_wt <- x$Cq[x$input == 1 & x$type == "+RT" & x$condition == "wildtype"]
        ref_rt_wt   <- x$Cq[x$input == 0 & x$type == "+RT" & x$condition == "wildtype"]
        adj_input_rt_wt <- input_rt_wt - 4.322                          # adjust input to 20 PCR cycles (log2=4.322)
        adj_input_rt_wt_minus_ref <- adj_input_rt_wt - mean(ref_rt_wt)  # substract reference RNA
        percent_input_rt_wt <- 100 * 2^adj_input_rt_wt_minus_ref        # calculate the percentage

        # calculate percentage input of -RT & wildtype condition
        #input_wt <- x$Cq[x$input == 1 & x$type == "-RT" & x$condition == "wildtype"]
        input_wt <- input_rt_wt   # Lara uses always +RT input!
        ref_wt   <- x$Cq[x$input == 0 & x$type == "-RT" & x$condition == "wildtype"]
        adj_input_wt <- input_wt - 4.322                       # adjust input to 20 PCR cycles (log2=4.322)
        adj_input_wt_minus_ref <- adj_input_wt - mean(ref_wt)  # substract reference RNA
        percent_input_wt <- 100 * 2^adj_input_wt_minus_ref     # calculate the percentage

        # calculate percentage input of +RT & TAP condition
        input_rt_tap <- x$Cq[x$input == 1 & x$type == "+RT" & x$condition == "TAP"]
        ref_rt_tap   <- x$Cq[x$input == 0 & x$type == "+RT" & x$condition == "TAP"]
        adj_input_rt_tap <- input_rt_tap - 4.322                          # adjust input to 20 PCR cycles (log2=4.322)
        adj_input_rt_tap_minus_ref <- adj_input_rt_tap - mean(ref_rt_tap) # substract reference RNA
        percent_input_rt_tap <- 100 * 2^adj_input_rt_tap_minus_ref        # calculate the percentage

        # calculate percentage input of -RT & TAP condition
        #input_tap <- x$Cq[x$input == 1 & x$type == "-RT" & x$condition == "TAP"]
        input_tap <- input_rt_tap   # Lara uses always +RT input!
        ref_tap   <- x$Cq[x$input == 0 & x$type == "-RT" & x$condition == "TAP"]
        adj_input_tap <- input_tap - 4.322                       # adjust input to 20 PCR cycles (log2=4.322)
        adj_input_tap_minus_ref <- adj_input_tap - mean(ref_tap) # substract reference RNA
        percent_input_tap <- 100 * 2^adj_input_tap_minus_ref     # calculate the percentage

        # do the barplot
        # sem is based on the sd formula as described in: 
        mean_wt  <- mean(percent_input_rt_wt - mean(percent_input_wt))
        sd_wt  <- sqrt(sd(percent_input_rt_wt)^2 + sd(percent_input_wt)^2)
        sem_wt <- sd_wt / sqrt(sum(length(percent_input_rt_wt) + length(percent_input_wt)))
        mean_tap <- mean(percent_input_rt_tap - mean(percent_input_tap))
        sd_tap  <- sqrt(sd(percent_input_rt_tap)^2 + sd(percent_input_tap)^2)
        sem_tap <- sd_tap / sqrt(sum(length(percent_input_rt_tap) + length(percent_input_tap)))

#        y   <- c(mean_wt, mean_tap)
#        sem <- c(sem_wt, sem_tap)
#        names(y) <- c("wild-type", paste0(x$binder[1], "-TAP"))
#        b <- barplot(y,
#                     ylim=c(min(c(0, y - sem)), max(y + sem)),
#                     ylab="(%input +RT) - (%input -RT)",
#                     main=x$loop[1])
#        arrows(b, y + c(sem_wt, sem_tap), b, y - c(sem_wt, sem_tap), angle=90, code=3, length=.1)
        df <- data.frame(mean=c(mean_wt, mean_tap),
                         sd  =c(sem_wt  , sem_tap  ),
                         cond=c("wildtype", paste0(x$binder[1], "-TAP")))
        df$cond <- factor(df$cond, levels=unique(df$cond))
        df.replicates <- data.frame(value=c(percent_input_rt_wt  - mean(percent_input_wt),
                                            percent_input_rt_tap - mean(percent_input_tap)),
                                    cond =c(rep("wildtype", length(percent_input_rt_wt)),
                                            rep(paste0(x$binder[1], "-TAP"), length(percent_input_rt_tap))))
        df.replicates$cond <- factor(df.replicates$cond, levels=unique(df.replicates$cond))

        p <- ggplot(df, aes(x=cond, y=mean)) +
                 geom_bar(stat="identity", color="black", size= 1) +
                 geom_point(data=df.replicates, mapping=aes(x=cond, y=value), alpha=.8, size=3) +
                 geom_errorbar(aes(ymin=mean, ymax=mean+sign(mean)*sd), width=0.25, size=1) +
                 ggtitle(x$loop[1]) + 
                 xlab("") + 
                 ylab("") +
                 theme_classic() +
                 theme(plot.title=element_text(hjust=0.5, size=22), 
                       axis.title.x = element_text(size=18),
                       axis.title.y = element_text(size=18),
                       text = element_text(size=20),
                       axis.text.x = element_text(angle=90, hjust=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

        p
    })
})

p <- unlist(p, recursive=FALSE)
p[[1]] <- p[[1]] + ylab("(%input +RT) - (%input -RT)")
pl <- gridExtra::marrangeGrob(grobs=p, nrow=1, ncol=4)
ggsave("rip.pdf", pl)

#dev.off()
