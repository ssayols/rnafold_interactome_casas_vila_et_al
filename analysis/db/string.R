x <- read.delim(gzfile("string.txt.gz"), sep=" ")
x$protein1 <- gsub("^4932\\.", "", x$protein1)
x$protein2 <- gsub("^4932\\.", "", x$protein2)
y <- reshape2::dcast(x, protein1 ~ protein2, value.var="combined_score")
rownames(y) <- y[,1]
y <- as.matrix(y[,-1])
write.csv(y, file="string.csv")
