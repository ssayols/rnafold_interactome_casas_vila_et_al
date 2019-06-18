##############################
##
## getseq: donat un fitxer de posicions "chr ini fi nom strand"
## obtenir les sequencies fasta de saccer2 i crear un fitxer multifasta
##
##############################
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
genome <- BSgenome.Scerevisiae.UCSC.sacCer2

feat <- read.delim("feat.csv",head=F)

for(i in 1:nrow(feat)) {

	print(i)

	# agafar la sequencia
	s <- getSeq(genome,feat[i,1],feat[i,2],feat[i,3],as.character=T)
	if(feat[i,5] == "-") {	# fer el RC si el gen esta a la negativa
		s <- chartr("ACTG","TGAC",paste(rev(unlist(strsplit(s,""))),collapse=""))
	}

	# partir en trossos de 80 caracters
	s <- sapply(seq(1,nchar(s),by=80),function(i) {
			substr(s,i,min(c(i+79,nchar(s))))
	})

	# escriure al fitxer multifasta de sortida
	cat(paste(">",paste(feat[i,4:6],collapse=" ")),file="feat.fasta",append=T,fill=T)
	sapply(s,cat,file="feat.fasta",append=T,fill=T)
}
