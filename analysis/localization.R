########################
##
## Localization analysis taking Weissman's GFP data
## Validation of candidates from other yeast RBPs from various authors
##
########################
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(grid)
library(UpSetR)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/localization.pdf", onefile=TRUE)

compartments.all <- c("ambiguous", "mitochondrion", "vacuole", "spindle.pole", "cell.periphery",
                      "punctate.composite", "vacuolar.membrane", "ER", "nuclear.periphery",
                      "endosome", "bud.neck", "microtubule", "Golgi", "late.Golgi",
                      "peroxisome", "actin", "nucleolus", "cytoplasm", "ER.to.Golgi",
                      "early.Golgi", "lipid.particle", "nucleus", "bud")
compartments <- c("nucleus", "nucleolus", "cytoplasm", "mitochondrion")
mitchell       <- read.delim("./analysis/db/yeastRBP_Mitchell.txt"      , comment="#")  # former Parker
matia_gonzalez <- read.delim("./analysis/db/yeastRBP_Matia_Gonzalez.txt", comment="#")  # former Gerber
shchepachev    <- read.delim("./analysis/db/yeastRBP_Shchepachev.txt"   , comment="#")  # former Rappsilber
beckmann       <- read.delim("./analysis/db/yeastRBP_Beckmann.txt"      , comment="#")
loc            <- read.delim("./analysis/db/localization.txt"           , comment="#")
binders        <- read.csv("./analysis/heatmap_1+1_threshold.csv")                      # our list of binders

# deduplicate binders
binders <- by(binders, 1:nrow(binders), function(x) {
    data.frame(X=unlist(strsplit(x$X, ";")),
               x[,!grepl("^(X|name)", colnames(x))],    # gene name is no more valid for XXX;YYY rows
               row.names=NULL)
})
binders <- do.call(rbind, binders)

# validation against Mitchell & Matia_Gonzalez's db
img <- venn.diagram(list(Mitchell=mitchell$Systematic.Name, Matia_Gonzalez=matia_gonzalez$ORF, Casas=binders$X),
                    fill=brewer.pal(3,"Set1"), alpha=.3, lwd=0, cex=1.5, cat.cex=1.5, filename=NULL)
grid.newpage()
grid.draw(img)

# validation against Mitchell & Matia_Gonzalez & Shchepachev's db
img <- venn.diagram(list(Mitchell=mitchell$Systematic.Name, Matia_Gonzalez=matia_gonzalez$ORF, Shchepachev=shchepachev$ORF, Casas=binders$X),
                    fill=brewer.pal(4,"Set1"), alpha=.3, lwd=0, cex=1.5, cat.cex=1.5, filename=NULL)
grid.newpage()
grid.draw(img)

# same, as upset
upset(fromList(list(Mitchell=mitchell$Systematic.Name, Matia_Gonzalez=matia_gonzalez$ORF, Shchepachev=shchepachev$ORF, Beckmann=beckmann$ENSEMBL.ID, Casas=binders$X)))

# barplot
df <- data.frame(compartment=compartments,
                 counts=sapply(compartments, function(x) sum(binders$X %in% loc$yORF[loc[, x]])))

p <- ggplot(df, aes(x=compartment, y=counts)) +
    geom_bar(stat="identity") +
    theme_bw()
print(p)

# venn
img <- venn.diagram(lapply(compartments, function(x) binders$X[binders$X %in% loc$yORF[loc[, x]]]),
                    category.names=compartments, fill=brewer.pal(4,"Set1"), alpha=.3, lwd=0, cex=1.5,
                    cat.cex=1.5, filename=NULL)
grid.newpage()
grid.draw(img)

# same, but as an Upset
x <- as.matrix(loc[loc$yORF %in% binders$X, compartments.all])
mode(x) <- "integer"
upset(as.data.frame(x))

# save output
binders$loc <- sapply(binders$X, function(x) { 
    paste(compartments[sapply(compartments, function(y) x %in% loc$yORF[loc[, y]])], collapse="; ")
})
write.csv(binders, file="./analysis/localization.csv", row.names=F)

dev.off()
