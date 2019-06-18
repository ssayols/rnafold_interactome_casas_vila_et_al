#############################
##
## Determine whether a construct is on an exon
##
## Got the annotation from SGD, release R62 as described in Weissman's paper
## https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R62-1-1_20090218.tgz
##
#############################
library(S4Vectors)  # zipup
library(rtracklayer)
library(GenomicRanges)

# read constructs and annotation
constructs <- makeGRangesFromDataFrame(read.csv("constructs.csv"), keep.extra.columns=TRUE)
annotation <- import.gff("saccharomyces_cerevisiae_R62-1-1_20090221.gff")
annotation <- annotation[! annotation$type %in% c("chromosome", "gene", "region")]

# check overlaps and add annotation to the constructs
i <- as.data.frame(findOverlaps(constructs, annotation))
i <- lapply(split(i, i$queryHits), function(x) {
    data.frame(row=x$queryHits[1],
               ann=paste(apply(as.data.frame(annotation[x$subjectHits]), 1, function(x) paste(x$type,
                                                                                              x$seqnames,
                                                                                              x$start,
                                                                                              x$end,
                                                                                              x$gene,
                                                                                              sep="-")),
                         collapse=" // "),
               ann2=paste(as.data.frame(annotation[x$subjectHits])$Note, collapse=" // ")
    )
})
i <- do.call(rbind, i)

# add annotation and save constructs table
constructs$ann <- NA
constructs$ann[i$row] <- i$ann
constructs$ann2 <- NA
constructs$ann2[i$row] <- i$ann2
write.csv(constructs, "constructs_on_exon.csv")
