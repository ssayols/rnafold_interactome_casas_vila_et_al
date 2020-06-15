###########################
##
## Plot dms signal of the constructs, ordered by binder
##
##   1-Read input data
##   2-Compute the in vivo/vitro data for each structure, and plot them
##   3-Plot in-vitro DMS signal of the constructs grouped by binder
##   4-Run the locarna structure+sequence aligner of ultraconserved sequences (CARNA)
##     with the structures grouped by binder
##   5-Take the results from the previous step, and replot(3) based on the CARNA realignments
##
###########################
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(parallel)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
library(Biostrings)

RATIO_THRES=1       # enrichment ratio used to detect binders
STICKY_THRES=10     # max constructs a protein should bind
CORES=4
GENOME=BSgenome.Scerevisiae.UCSC.sacCer2
LOCARNA="/fsimb/groups/imb-buttergr/=EVOREG/rnalib/bin/locarna/bin/mlocarna"
CARNA="/fsimb/groups/imb-buttergr/=EVOREG/rnalib/bin/carna/bin/carna"
RUN_LOCARNA=FALSE

##
## 1-read input data
##
# read and parse DMS data. Format is wiggle
dms.vivo.minus  <- import.wig("./DMS-seq data/Silvi/May2013_VivoAllextra_minus.txt")
dms.vivo.plus   <- import.wig("./DMS-seq data/Silvi/May2013_VivoAllextra_plus.txt")
dms.vitro.minus <- import.wig("./DMS-seq data/Silvi/May2013_VitroAllextra_minus.txt")
dms.vitro.plus  <- import.wig("./DMS-seq data/Silvi/May2013_VitroAllextra_plus.txt")
dms.denat.minus <- import.wig("./DMS-seq data/Silvi/May2013_DenaturedAllextra_minus.txt")
dms.denat.plus  <- import.wig("./DMS-seq data/Silvi/May2013_DenaturedAllextra_plus.txt")

# read input structures
structs <- read.csv("./constructs/constructs.csv")
structs <- with(structs, GRanges(seqnames=chr, ranges=IRanges(start, end), strand=strand, name=name))
structs$seq <- getSeq(GENOME, structs, as.character=T)

# read binders
binders <- read.csv("./analysis/heatmap_1+1_threshold.csv", row.names=1)

##
## 2-Compute the in vivo/vitro data for each structure, and plot them
##
get.dms.data <- function(x) {
    # get in vivo/vitro data
    dms.vivo  <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.vivo.minus, x),
                                      "+"=subsetByOverlaps(dms.vivo.plus , x)))
    dms.vitro <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.vitro.minus, x),
                                      "+"=subsetByOverlaps(dms.vitro.plus , x)))
    dms.denat <- as.data.frame(switch(as.character(strand(x)),
                                      "-"=subsetByOverlaps(dms.denat.minus, x),
                                      "+"=subsetByOverlaps(dms.denat.plus , x)))

    # normalize to the most reactive base
    dms.vivo$nscore  <- dms.vivo$score  / max(dms.vivo$score)
    dms.vitro$nscore <- dms.vitro$score / max(dms.vitro$score)
    dms.denat$nscore <- dms.denat$score / max(dms.denat$score)

    # merge by pos
    dms <- Reduce(function(x, y) merge(x, y, by=c("start", "end"), all=T), list(dms.vivo, dms.vitro, dms.denat))
    dms <- dms[, c("seqnames", "start", "nscore.x", "nscore.y", "nscore")]
    colnames(dms) <- c("chr", "pos", "vivo", "vitro", "denat")

    # plot the DMS signal along the structure
    x <- melt(dms, id.vars=c("chr", "pos"))
    mid <- min(x$pos) + round((max(x$pos) - min(x$pos)) / 2)

    p <- ggplot(x, aes(x=pos, y=value, fill=variable)) +
            geom_bar(position="dodge", stat="identity") +
            scale_fill_discrete("", guide=FALSE) +
            scale_x_continuous(x$chr[1], breaks=c(min(x$pos), mid, max(x$pos))) +
            ylab("Normalized DMS signal") +
            theme_minimal() +
            facet_wrap(~variable, dir="v")

    # return
    list(p, dms)
}

# and plot
dms.plots <- mclapply(structs, get.dms.data, mc.cores=CORES)
names(dms.plots) <- structs$name

cairo_pdf("./analysis/DMSbinders_constructs.pdf", onefile=TRUE)
mapply(function(x, y) { print(x[[1]] + ggtitle(y)); invisible(0) }, dms.plots, names(dms.plots))
dev.off()

##
##   3-Plot in-vitro DMS signal of the constructs grouped by binder
##
sel <- apply(binders[, -1], 1, function(x){
    y <- sum(x > RATIO_THRES)
    y > 0 & y <= STICKY_THRES
})

dms.binders <- apply(binders[sel, ], 1, function(x) {

    # select the DMS from the constructs where the protein is enriched
    gene  <- x[1]
    cons  <- names(x[-1])
    ratio <- as.numeric(x[-1])
    dms.sel <- dms.plots[cons[ratio > RATIO_THRES]]
    dms.sel <- lapply(dms.sel, function(x) melt(x[[2]], id.vars=c("chr", "pos")))

    # plot the DMS signal along the structure
    df <- melt(dms.sel, measure.vars="value")
    df <- df[df$variable == "vitro",]

    p <- ggplot(df, aes(x=pos, y=value, fill=L1)) +
            geom_bar(position="dodge", stat="identity") +
            scale_fill_discrete("", guide=FALSE) +
            ylab("Normalized DMS signal") +
            theme_minimal() +
            ggtitle(gene) +
            facet_wrap(~L1, ncol=1, dir="v", scales="free")
   
    # and return everything
    list(gene=gene, cons=cons, ratio=ratio, dms.sel=dms.sel, df=df, ggplot=p)
})

cairo_pdf("./analysis/DMSbinders.pdf", width=8.267, height=11.692, onefile=TRUE)
x <- lapply(dms.binders, function(x) print(x$ggplot))
dev.off()

##
##   4-Run the locarna structure+sequence aligner of ultraconserved sequences (CARNA)
##     with the structures grouped by binder
##
if(!dir.exists("./analysis/DMSbinders_msa")) dir.create("./analysis/DMSbinders_msa")
x <- lapply(dms.binders, function(x) {

    # save in fasta format the structures for this binder
    x$gene <- make.names(x$gene)
    out  <- paste0("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/DMSbinders_msa/", x$gene, ".fa")  # output fasta file
    s    <- as.data.frame(structs[structs$name %in% x$cons[x$ratio > RATIO_THRES]])  # structures it binds to
    apply(s, 1, function(s) cat(">", s["name"], "\n", s["seq"], "\n", sep="", fill=FALSE, append=TRUE, file=out))   # >header\nseq\n

    # run the locarna structure+sequence aligner of ultraconserved sequences
    threads <- min(CORES, sum(x$ratio > RATIO_THRES))
    if(RUN_LOCARNA && threads > 1 && !file.exists(paste0(out, ".out"))) {
        cat("submitting LOCARNA for", basename(out), fill=T)
        cmd <- paste0("echo '", LOCARNA, " --threads=", threads, " --pw-aligner ", CARNA, " ", out, "' | ",
                      "bsub -J ", x$gene, " -W5:00 -n", threads, " -R'span[ptile=", threads, "]' -app Reserve1G -o ", out, ".out")
        system(cmd)
    } else {
        cat("skipping LOCARNA for", basename(out), fill=T)
    }
})

## WAIT HERE UNTIL ALL BSUB JOBS ARE FINISHED ##

##
##   5-Take the results from the previous step, and replot(3) based on the CARNA realignments
##
dms.binders <- lapply(dms.binders, function(x) {

    # if the binders just binds to 1 struct, plot the same without realignment (no locarna run)
    out <- paste0("/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/DMSbinders_msa/", x$gene, ".out/results/result.aln")
    if(sum(x$ratio > RATIO_THRES) <= 1 | !file.exists(out)) {
        x$ggplot_aln <- x$ggplot
    } else {
        # read the output alignment from locarna
        x$gene <- make.names(x$gene)
        aln  <- as.character(readRNAMultipleAlignment(out))

        # get the new coordinates of the structure after alignment
        newpos <- mapply(function(name, s) {
            ini <- start(structs[structs$name == name])
            x  <- data.frame(seq=unlist(strsplit(s, "")),
                             pos=0,
                             newpos=1:nchar(s))

            # add the alignment offset
            i     <- gregexpr("[AUCG]", s)[[1]]
            x$pos[i] <- ini:(ini + length(i) - 1)

            x
        }, names(aln), aln, SIMPLIFY=FALSE)

        # reposition the signal based on the anaysis
        df <- by(x$df, x$df$L1, function(x) {
            x$newpos <- newpos[[x$L1[1]]]$newpos[match(x$pos, newpos[[x$L1[1]]]$pos)]
            x
        })
        x$df <- do.call(rbind, df)

        # and merge all structures in a consensus struct merged by max DMS signal
        consensus.max <- by(x$df$value, x$df$newpos, mean, na.rm=T) # na.rm removes NA signals for in-vitro
        x$df <- rbind(x$df, data.frame(chr=x$df$chr[1],
                                       pos=as.integer(names(consensus.max)),
                                       variable="vitro",
                                       value=as.numeric(consensus.max) / max(consensus.max, na.rm=T), # normalize again
                                       L1="Consensus (normalized to the max average signal)",
                                       newpos=as.integer(names(consensus.max))))

        # change the ggplot data with the new coords
        x$ggplot_aln <- x$ggplot %+% x$df + aes(x=newpos, y=value, fill=L1) + facet_wrap(~L1, ncol=1, dir="v")
    }
    
    x
})

cairo_pdf("./analysis/DMSbinders(aligned).pdf", width=8.267, height=11.692, onefile=TRUE)
x <- lapply(dms.binders, function(x) print(x$ggplot_aln))
dev.off()
