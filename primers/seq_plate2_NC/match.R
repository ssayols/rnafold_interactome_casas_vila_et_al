##############################
##
## Buscar si una sequencia query esta dins les sequencies target
##
##############################
library("ShortRead")
library("rmarkdown")
setwd("~/IMB/nuria/seq_plate2_NC")
QUERY  <- "../feat.fasta"	# fitxer d'entrada amb les sequencies query

## llegir fitxer multifasta
wells  <- read.csv("yeastID_wellposition.csv")
query  <- readFasta(QUERY)	# a vector of objects
target <- lapply(list.files(pattern="77BB97.*\\.fas"),readFasta)	# a LIST of objects

# calcular reverse complimentari
RC <- function(s) {
	chartr("ACTG","TGAC",paste(rev(unlist(strsplit(s,""))),collapse=""))
}

# progress bar
progress.old <<- 0
progress <- function(current){
	if(current > progress.old) { # just print every 1% increase, not for every step
		cat(paste(c("\r[",rep("=",current),rep(" ",100-current),"]",current,"%"),collapse=""))
		flush.console()
		progress.old <<- current
	}
}

# taula de sortida
query_id <- gsub("^ ","",as.character(id(query)))
query_id.short <- gsub("^.+ [-+] (.+)$","\\1",query_id)
target_id <- unlist(lapply(target,function(x) gsub(".+- ID: (.+)- on.+","\\1",as.character(id(x)))))

df <- data.frame(oligoname=query_id,
				 oligoname.short=query_id.short,
				 gatc="",
				 filler="",
				 plate=wells$Plate[match(query_id.short,wells$Oligoname)],
				 position=wells$Position[match(query_id.short,wells$Oligoname)],
				 match="")

## loop per cada sequencia
for(i in 1:length(query)) {

	progress(round(100 * i / length(query)))
	
	# buscar tots els matchs	
	j <- lapply(target,function(target,query) {
		query  <- toupper(as.character(sread(query)))
		target <- toupper(as.character(sread(target)))
		forward <- regexpr(query,target)
		reverse <- regexpr(RC(query),target)
		return(ifelse(forward != -1,forward,
			   ifelse(reverse != -1,-1 * reverse,0)))
	},query[i])
	j <- unlist(j)
	
	# omplir la taula
	if(any(j != 0)) {

		pos_query  <- which(df$oligoname == gsub("^ ","",as.character(id(query[i]))))
		pos_target <- which(j != 0)

		# control de targets que tenen la mateixa sequencia (enviats 2+ cops a sequenciar a gatc)
		if(length(pos_target > 1)) {
			if(any(j[pos_target] > 0)) {
				pos_target <- pos_target[j[pos_target] > 0][1]	# reportem nomes el primer positiu
			}
			else {
				pos_target <- pos_target[1]	# reportem nomes el primer negatiu
			}
		}
		df[pos_query,"gatc"]  <- target_id[pos_target]
		df[pos_query,"match"] <- j[pos_target]

		# escriure MD i render
		if(j[pos_target] > 0) {
			target_seq <- as.character(sread(target[[pos_target]]))
			query_seq  <- as.character(sread(query[i]))
			out <- c(paste0("#### ",target_id[pos_target]," - ",df$oligoname[pos_query]," ####"),
					 paste0(substr(target_seq,1,j[pos_target] - 1),
							"**",query_seq,"**",
							substr(target_seq,j[pos_target] + nchar(query_seq),nchar(target_seq))))
			f <- file("tmp.Rmd")
			writeLines(out,f)
			close(f)
			render("tmp.Rmd","word_document",paste0(target_id[pos_target],".docx"))
			file.remove("tmp.Rmd")
		}
	}
}

write.csv(df,file="match.csv",row.names=F)
