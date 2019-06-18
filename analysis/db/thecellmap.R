library(parallel)

# read in the tables (ExE, NxE and NxN)
cl <- makeCluster(3)
f  <- c("SGA_ExN_NxE_clustered.txt", "SGA_ExE_clustered.txt", "SGA_NxN_clustered.txt")
cellmap <- parSapply(cl, f, function(f) {
    x <- read.delim(f)
    do.call(rbind, by(x, x[, 1], function(x) apply(x[, -1], 2, mean, rn.rm=T)))
})
stopCluster(cl)

##
## merge the 3 tables together. Assumptions:
##   -all cols in x[[2]] and x[[3]] are also in x[[1]], but 99.9999% NA in x[[1]]
##   -some rows in x[[2]] and x[[3] are also in x[[1]], and need to be added (merged)
##
merge.cellmap <- function(x, y) {
    # ExN_NxE (cellmap[[1]]) contains incomplete data about ExE and NxN. Add it
    common.cols <- colnames(x)[colnames(x) %in% colnames(y)]
    common.rows <- rownames(x)[rownames(x) %in% rownames(y)]
    x[common.rows, common.cols] <- y[common.rows, common.cols]

    # merge the non common rows
    y.unique <- y[!(rownames(y) %in% rownames(x)), common.cols]
    m <- matrix(rep(NA, nrow(y.unique) * (ncol(x) - length(common.cols))), nrow=nrow(y.unique))
    colnames(m) <- colnames(x)[!(colnames(x) %in% common.cols)]
    y.unique <- cbind(y.unique, m)
    x <- rbind(x, y.unique[, colnames(x)])

    x
}
cellmap.merged <- Reduce(merge.cellmap, cellmap)

write.csv(cellmap.merged, file="thecellmap.csv")
