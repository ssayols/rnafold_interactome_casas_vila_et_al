options(stringsAsFactors=FALSE)
library(seqinr)

setwd("/Volumes/groups/imb-buttergr/=EVOREG/rnalib/constructs/calculate frames CDS/")

# sequences and constructs
translation_table <- read.csv("../../constructs/constructs.csv")
translation_table$fasta_name <- paste(">",
                                      translation_table$name1,
                                      translation_table$strand,
                                      translation_table$name2)

feat <- seqinr::read.fasta("../../primers/feat.fasta")
p_fwd <- list("+1"=s2c(c("GTACGCCTCGAGCCGGACTCTAGAGGATCGAACCCTT"  )),     #2683 len 31 + 2 from construct's sequence
              "+2"=s2c(c("GTACGCCTCGAGCCCCGGACTCTAGAGGATCGAACCCTT")),     #2717 len 32 + 1 from construct's sequence
              "+0"=s2c(c("GTACGCCTCGAGCCCGGACTCAAGAGGATCGAACCCTT")))      #2718 len 33 + 0
p_rev <- list("+2"=s2c(c("AAGGGTTCGATCCCTACCGGTTACCTTAAGGTACG"  )),   #2612 len 35 = 33+2 = p_fwd+construct+2+33
              "+0"=s2c(c("AAGGGTTCGATCCCTACCGGTTACCCTTAAGGTACG" )),   #2613 len 36 = 36   = p_fwd+construct+0+36
              "+1"=s2c(c("AAGGGTTCGATCCCTACCGGTTACCCCTTAAGGTACG")))   #2614 len 37 = 36+1 = p_fwd+construct+1+36

# process sequences
out <- sapply(seq_along(feat), function(i) {
  
  name <- getAnnot(feat[[i]])
  aa   <- lapply(1:3, function(shift) getTrans(c(p_fwd[[shift]], getSequence(feat[[i]]))))
  
  # p_fwd controls the frame shift
  valid_shifts <- which(grepl("\\*", aa) == 0)
  if(length(valid_shifts) > 0) {
    
    # p_fwd + construct_sequence
    p1 <- p_fwd[[valid_shifts[1]]]
    fwd_seq <- c(p1, getSequence(feat[[i]]))
    len  <- getLength(fwd_seq)
    
    # p_rev controls the sequence length, to be multiple of 3
    p2 <- switch(as.character(len %% 3),
                 "0" = p_rev[["+0"]],
                 "1" = p_rev[["+1"]],
                 "2" = p_rev[["+2"]])
    fwd_seq_rev <- c(fwd_seq, p2)
  } else {
    p1 <- NA
    p2 <- NA
    fwd_seq_rev <- NA
  }
  
  # output
  data.frame(name   =translation_table$name[match(name, translation_table$fasta_name)],
             name2  =sub("^> ", "", name),
             gene   =translation_table$external_gene_name[match(name, translation_table$fasta_name)],
             pos    =translation_table$pos[match(name, translation_table$fasta_name)],
             length =getLength(getSequence(feat[[i]])),
             length2=getLength(fwd_seq_rev),
             frames =paste0(names(p_fwd)[valid_shifts], collapse="; "),
             p1     =paste0(p1, collapse=""),
             p2     =paste0(p2, collapse=""),
             seq    =paste0(getSequence(feat[[i]]), collapse=""),
             seq2   =paste0(fwd_seq_rev, collapse=""),
             seqAA  =paste0(getTrans(fwd_seq_rev), collapse="")
  )
})

write.csv(t(out), file="calculate_frames.2.csv", row.names = FALSE)
