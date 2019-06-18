##############################
##
## Fer primers a partir d'un fitxer multifasta
##
## 1-llegir fitxer multifasta
## 2-loop per cada sequencia:
##   2.1-agafar n=MINNUC del costat 5'. Agafar n=n+1 fins que TM(seq_n) > TM
##   2.2-agafar n=MINNUC del costat 3'. Agafar n=n+1 fins que TM(RC(seq_n)) > TM
##   2.3-output: nom;primer_5';TM_5';primer_3';TM_3'
##
## TM=4*(C+G)+2*(A+T)
##
##############################
library("ShortRead")
FASTA <- "feat.fasta"	# fitxer d'entrada amb les sequencies dels gens
MINL  <- 15
TM    <- 58
OUT   <- "doprim.csv"	# fitxer de sortida amb els primers

## 1-llegir fitxer multifasta
fasta <- readFasta(FASTA)

# calcul de la temperatura de fusio
tm <- function(s) {

	nC <- length(gregexpr("C",s)[[1]])
	nG <- length(gregexpr("G",s)[[1]])
	nA <- length(gregexpr("A",s)[[1]])
	nT <- length(gregexpr("T",s)[[1]])

	return(4*(nC+nG)+2*(nA+nT))
}

# calcular reverse complimentari
RC <- function(s) {
	chartr("ACTG","TGAC",paste(rev(unlist(strsplit(s,""))),collapse=""))
}

## 2-loop per cada sequencia:
for(i in 1:length(fasta)) {

	s <- as.character(sread(fasta[i]))	# sequencia
	x <- as.character(id(fasta[i]))		# id

	##   2.1-agafar n=MINNUC del costat 5'. Agafar n=n+1 fins que TM(seq_n) > TM
	n  <- MINL
	t5 <- 0
	while(t5 < TM) {
		s5 <- substr(s,1,n)
		t5 <- tm(s5)
		n  <- n + 1
	}

	##   2.2-agafar n=MINNUC del costat 3'. Agafar n=n+1 fins que TM(RC(seq_n)) > TM
	n  <- MINL
	t3 <- 0
	while(t3 < TM) {
		s3 <- RC(substr(s,nchar(s) - n,nchar(s)))	# idem, pero Reverse Complimentari
		t3 <- tm(s3)
		n  <- n + 1
	}

	##   2.3-output: nom;primer_5';TM_5';primer_3';TM_3'
	cat(paste(x,s5,t5,nchar(s5),s3,t3,nchar(s3),s,sep=";"),file=OUT,fill=T,append=T)
}
