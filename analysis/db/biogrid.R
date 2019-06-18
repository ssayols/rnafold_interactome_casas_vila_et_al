x <- read.delim(gzfile("biogrid.txt.gz"), sep="\t")
x$Score <- as.numeric(ifelse(x$Score == "-", 1, x$Score))   # give some value to NA scores
x <- x[, c("Systematic.Name.Interactor.A", "Systematic.Name.Interactor.B", "Score")]
x <- x[x$Systematic.Name.Interactor.A != "-" & x$Systematic.Name.Interactor.B != "-", ]
x <- unique(x, MARGIN=1)    # remove duplicated interactions from different studies, and get the best score
x <- by(x, paste(x$Systematic.Name.Interactor.A, x$Systematic.Name.Interactor.B), function(x) x[which.max(abs(x$Score)), ])
x <- data.table::rbindlist(x)
y <- reshape2::dcast(x, Systematic.Name.Interactor.A ~ Systematic.Name.Interactor.B, value.var="Score")
rownames(y) <- y[,1]
y <- as.matrix(y[,-1])
write.csv(y, file="biogrid.csv")
