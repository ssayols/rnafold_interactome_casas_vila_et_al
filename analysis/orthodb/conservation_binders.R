##############################
##
## Check conservation of binders in other species, namely: 
##
##   * Taxonomy level 4893: Saccharomycetaceae
##   * Species ids:
##     - 1064592: Naumovozyma castellii
##     - 1160507: Saccharomyces arboricola
##     - 226230: Saccharomyces kudriavzevii
##     - 559292: Saccharomyces cerevisiae S288c
## 
## All work is done on the data available on Orthodb v.9, and we're heavily 
## restricted with the species to use for our comparison (though these 3 are ok
## for a not very serious conservation analysis).
##
## Steps:
##   1-Recreate the Orthodb v.9 database to get the protein sequence of the ortholog proteins in other species
##   2-Perform MSA (ClustalO, using its implementation in the Bioconductor's msa packag)
##     for each orthology group (which in principle represent 1 SacCer protein)
##   3-get some sort of metric and plot against a background set of non-binders
##
####################################
library(DBI)
library(biomaRt)
library(Biostrings)
library(msa)
library(ggplot2)
data(BLOSUM62)

set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")
cairo_pdf("./analysis/orthodb/conservation_binders.pdf", onefile=TRUE)

SPECIES <- data.frame(ncbi_tax_id=c(1064592,
                                    1160507,
                                    226230,
                                    559292),
                      scientific_name=c("Naumovozyma castellii",
                                        "Saccharomyces arboricola",
                                        "Saccharomyces kudriavzevii",
                                        "Saccharomyces cerevisiae"))

ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "uniprot_swissprot"),
             mart=useMart("ENSEMBL_MART_ENSEMBL", host="oct2016.archive.ensembl.org", dataset="scerevisiae_gene_ensembl"))

db <- dbConnect(RSQLite::SQLite(), dbname="./analysis/orthodb/ortho.db")

##
## get all orthologs in other species and do msa
##
binders <- unlist(strsplit(read.csv("./analysis/heatmap_1+1_threshold.csv")$X, ";"))
binders <- data.frame(gene_id   =binders,
                      gene_name =ann$external_gene_name[match(binders, ann$ensembl_gene_id)],
                      uniprot_id=ann$uniprot_swissprot[match(binders, ann$ensembl_gene_id)])
binders$gene_name <- ifelse(binders$gene_name == "", binders$gene_id, binders$gene_name)

scores <-
    lapply(list(binders=binders$uniprot_id,
                bg=sample(ann$uniprot_swissprot, 1000)),#length(binders))),
           function(x) {
    
        query <- paste(sep="\n",
            # Get sequences of orthologous genes for the 4 species up here
            "SELECT a.organism_tax_id, a.uniprot_id, a.description, b.OG_unique_id, c.seq",
            "  FROM genes AS a",
            "  JOIN OG2genes AS b ON a.Ortho_DB_gene_id=b.Ortho_DB_gene_id",
            "  JOIN fasta_fungi AS c ON a.uniprot_id=c.uniprot_id",
            " WHERE b.OG_unique_id IN (",
            #    Get orthologous groups of our gene of interest, calculated at the Saccharomycetaceae level
            "   SELECT b.OG_unique_id",
            "     FROM genes AS a",
            "     JOIN OG2genes AS b ON a.Ortho_DB_gene_id=b.Ortho_DB_gene_id",
            "     JOIN OGs AS c ON b.OG_unique_id=c.OG_unique_id",
            "    WHERE a.uniprot_id IN (",
            paste(paste0("'", sort(unique(x)), "'"), collapse=","),
            "    )",
            "      AND c.level_tax_id_on_which_the_group_was_built=4893",
            " )",
            "   AND a.organism_tax_id IN (1064592, 226230, 559292, 1160507)",
            " ORDER BY b.OG_unique_id, a.organism_tax_id, a.uniprot_id;")
    
        x <- dbGetQuery(db, query)
    
        ##
        ## MSA each Orthology Group together, and return the mean conservation score of each residue
        ##
        tryCatch({
            sapply(split(x, x$OG_unique_id), function(x) {   # dont run in parallel, it messes up with temp files!
                if(length(unique(x$organism_tax_id)) < 2) return(NA)

                # annotate what will be the fasta headers
                # > <uniprot_id> <organism> <ensembl_gene_id of the SacCer binder>
                # > G0V636 Naumovozyma castellii PAB1
                binder <- paste(ann[match(x$uniprot_id[x$organism_tax_id == "559292"], ann$uniprot_swissprot), ], collapse="_")
                x$fasta_head <- paste(x$uniprot_id,
                                      SPECIES$scientific_name[match(x$organism_tax_id, SPECIES$ncbi_tax_id)], 
                                      binder)
        
                # build the AAStringSet with all sequences
                seqs <- x$seq
                names(seqs) <- x$fasta_head
                seqs <- AAStringSet(seqs)
        
                # do msa (ClustalW, all defaults [substitution matrix, etc.])
                m <- msa(seqs)
                log(mean(msaConservationScore(m, BLOSUM62), na.rm=TRUE))
            })
        }, error=function(e) NA)
    })


##
## Some plotting
##

do_plot <- function(scores, main="") {

    # as a barplot (mean + error bars)
    df <- data.frame(condition=names(scores),
                     avg      =sapply(scores, mean, na.rm=TRUE),
                     dev      =sapply(scores, sd  , na.rm=TRUE))

    p <- ggplot(df, aes(x=condition, y=avg, fill=condition)) +
             geom_bar(stat="identity", position="dodge") +
             geom_errorbar(aes(ymin=avg-dev, ymax=avg+dev), width=.2, position=position_dodge(.9)) +
             scale_fill_manual("", values=c("#E69F00", "#56B4E9")) +
             ylab("mean conservation score (log)") + xlab("") +
             ggtitle(main) +
             theme_bw()
    print(p)
#    readline("continue...")

    # as a boxplot
    x <- reshape::melt(scores)
    p <- ggplot(x, aes(x=L1, y=value, fill=L1)) +
             geom_violin(aes(color=L1), trim=FALSE, show.legend=FALSE) +
             geom_boxplot(width=.1) +
             scale_fill_manual("", values=c("#E69F00", "#56B4E9")) +
             scale_color_manual("", values=c("#E69F00", "#56B4E9")) +
             ylab("mean conservation score (log)") + xlab("") +
             ggtitle(main) +
             theme_bw()
    p <- ggpval::add_pval(p, pairs=list(c(1, 2)))
    print(p)
#    readline("continue...")

    # density distributions
    pal <- c("red", "blue")
    x <- lapply(scores, density, na.rm=TRUE)
    all.x <- do.call(c, lapply(x, function(x) x$x))
    all.y <- do.call(c, lapply(x, function(x) x$y))
    plot(0, xlim=c(min(all.x), max(all.x)), ylim=c(min(all.y), max(all.y)), xlab="", ylab="", type="n", bty="n", main=main)
    Map(x, col=pal, f=lines)
    legend("topright", legend=names(x), fill=pal, bty="n")
#    readline("continue...")

    invisible(0)
}

# conservation scores considering the proteins that don't have an ortholog at all (conservation score = 0)
do_plot(scores,
        "*not* considering proteins that don't have orthologs in other species")
do_plot(lapply(scores, function(x) ifelse(is.na(x), 0, x)),
        "considering proteins that don't have orthologs in other species\n(conservation score = 0)")

dev.off()
