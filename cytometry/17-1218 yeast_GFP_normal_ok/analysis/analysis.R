setwd("/Volumes/Macintosh HD/Users/ncasasvi/Documents/PROJECTS/_mRNA fold quantitative proteomics/FUNCTIONAL FOLLOW-UP/functional UTR screen/CYTOMETRY/Results/17-1218 yeast_GFP/analysis/")
library(ggplot2)
library(WriteXLS)
results <- read.delim("results.txt", stringsAsFactors = FALSE)

pdf("results.pdf") #width=7, height=6)

res <- lapply(split(results, results$RBP), function(results.RBP) {

  results.RBP$fold <- factor(results.RBP$fold,
                             levels=c(unique(results.RBP$fold[results.RBP$type == "experiment"]),
                                      unique(results.RBP$fold[results.RBP$type != "experiment"])))
  
  p <- ggplot(results.RBP, aes(x=fold, y=mean, fill=type)) +
    geom_bar(stat="identity", position=position_dodge(), color="black", size= 1) +
    scale_fill_manual(values=c("#999999", "#E69F00")) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), position=position_dodge(), width = 0.25, size = 1) +
    ggtitle(results.RBP$RBP[1]) + 
    xlab("RNA fold") + 
    ylab("GFP ko/GFP wildtype") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, size=22), 
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          text = element_text(size=20),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  print(p)
  
  x <- reshape::melt.data.frame(results.RBP[, 2:5]) # columnes fold + rep[123]
  summary(a <- aov(x$value ~ x$fold))
  as.data.frame(TukeyHSD(a)$`x$fold`)
})

WriteXLS(res, ExcelFileName = "anova.xls", row.names = TRUE)

dev.off()
