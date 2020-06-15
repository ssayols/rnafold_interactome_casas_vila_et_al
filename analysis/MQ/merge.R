#############################
##
## filter and merge PG from the 2 plates
##
#############################
options(stringsAsFactors=F)
library(rMQanalysis)
library(parallel)

setwd("/fsimb/groups/imb-buttergr/=EVOREG/rnalib")

PG <- list("./analysis/MQ/proteinGroups_plate1.txt",
           "./analysis/MQ/proteinGroups_plate2.txt")

pg <- mclapply(PG, function(x) {
    filterIdentifiedProteins(filterWholeDataset(read.delim(x)))
}, mc.cores=2)

pg <- merge(pg[[1]], pg[[2]], by="Majority.protein.IDs", all=T)

write.csv(pg, file="/fsimb/groups/imb-buttergr/=EVOREG/rnalib/analysis/MQ/proteinGroups.txt", row.names=F)
