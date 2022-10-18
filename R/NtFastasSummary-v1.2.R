
library(Biostrings)
library(stringr)
library(data.table)

###  Read amplicon aligned sequences 
read.ampl.seqs <- function(flnm,mnr=1)
{
  seqs <- as.character(readDNAStringSet(flnm))
  IDstr <- names(seqs)
  n <- length(IDstr)
  nms <- character(length=n)
  nr <- integer(n)
  for(j in 1:n)
  { strs <- strsplit(IDstr[j],split="\\|")[[1]]
    nms[j] <- strs[1]
    if(length(strs)>1)
	  nr[j] <- as.numeric(strs[2])
  }
  IDs <- data.frame(ID=nms,nseqs=nr,stringsAsFactors=FALSE)
  nall <- nrow(IDs)
  tnr <- sum(IDs$nseqs)
  return(list(IDs=IDs,seqs=seqs,nall=nall,tnr=tnr))
}


SortByMutations <- function(bseqs,nr)
{
  master <- bseqs[which.max(nr)]
  ##  Determinar diferencias respecte a la máster
  bseqs <- DNAStringSet(bseqs)
  psa <- pairwiseAlignment(pattern=bseqs,subject=master)
  nm <- nmismatch(psa)
  tnm <- table(nm)
  ##  Ordenar por número de mutaciones
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]
  ##  Número de orden dentro de cada número de mutaciones
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  ##  Ordenar por frecuencia descendiente dentro de cada número de mutaciones
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }
  ##  Calcular frecuencia relativa
  frq <- round(nr/sum(nr)*100,2)
  ## Nombre completo para cada haplotipo
  nms <- paste("Hpl",nm,sprintf("%04d",isq),sep="_")
  ##  Cabecera fasta con nombre, número de reads y frecuencia relativa
  names(bseqs) <- nms
  list(bseqs=as.character(bseqs),nr=nr,nm=nm)
}


###  Read amplicon aligned sequences with abundance data, filter
###    at a minimum given abundance, and sort by mutations and
###    abundance.
get.data <- function(flnm,min.pct=0.1)
{ lst <- read.ampl.seqs(flnm)
  fl <- lst$IDs$nseqs/sum(lst$IDs$nseqs)*100 >= min.pct
  seqs <- lst$seqs[fl]
  names(seqs) <- lst$IDs$ID[fl]
  nr <- lst$IDs$nseqs[fl]
  lst <- SortByMutations(seqs,nr)
  list(seqs=lst$bseqs,nr=lst$nr,nm=lst$nm)
}

#------------------------------------------------------------------------#

primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
rownames(primers) <- primers$Ampl.Nm

fl.patt <- "\\.fna$"

flnms <- list.files(path=ntDir,patt=fl.patt)
parts <- t(sapply(flnms,function(str) strsplit(str,split="\\.")[[1]]))
anms <- parts[,4]
snms <- parts[,1]

n <- length(flnms)
rprt <- data.frame(ID=snms,Ampl=anms,Hpl=integer(n),MaxMut=integer(n),
                   Master=numeric(n),Reads=numeric(n),stringsAsFactors=FALSE)

for(i in 1:n)
{
   lst <- get.data(file.path(ntDir,flnms[i]),0)   
   rprt$Hpl[i] <- length(lst$nr)
   rprt$MaxMut[i] <- max(lst$nm)
   rprt$Master[i] <- round(max(lst$nr)/sum(lst$nr)*100,2)
   rprt$Reads[i] <- sum(lst$nr)
}

o <- order(rprt$Ampl,rprt$ID)

sink(file.path(repDir,"NtFastasSummary.txt"))
cat("\nFasta files summary\n\n")
print(rprt[o,])
cat("\nFinal coverage, global summary\n\n")
print(summary(rprt$Reads))
cat("\nFinal coverage, by amplicon summary\n\n")
print(tapply(rprt$Reads,rprt$Ampl,summary))
sink()

library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")

pdf(file.path(repDir,"FinalCovBoxplots.pdf"),paper="a4",width=5,height=10)
par(mfrow=c(2,1))

y <- log10(rprt$Reads)
boxplot(y~rprt$Ampl,border="gray",ylab="Log10 #reads",outline=FALSE,
        ylim=range(y),las=2)
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),y,pch="+",cex=0.8)	
title(main="Final coverage by amplicon")	
abline(h=log10(10000),col="red",lty=4)

o <- order(rprt$Reads)
y <- log10(rprt$Reads[o])
g <- factor(rprt$Ampl[o])
icol <- pal[as.integer(g)]
plot(1:length(y),y,type="h",col=icol,ylab="log10 #reads",xlab="sorted samples")
legend("topleft",lwd=2,col=pal,legend=levels(g),cex=0.8)
abline(h=log10(10000),col="red",lty=4)

boxplot(rprt$Reads~rprt$Ampl,border="gray",ylab="#reads",outline=FALSE,
        ylim=range(rprt$Reads))
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),rprt$Reads,pch="+",cex=0.8)		
abline(h=10000,col="red",lty=4)

y <- log10(rprt$Reads)
boxplot(y~rprt$Ampl,border="gray",ylab="Log10 #reads",outline=FALSE,
        ylim=range(y),las=2)
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),y,pch="+",cex=0.8)	
title(main="Final coverage by amplicon")	
abline(h=log10(10000),col="red",lty=4)

dev.off()

file.copy(file.path(repDir,"FinalCovBoxplots.pdf"),
          file.path(resultsDir,"FinalCovBoxplots.pdf"),
		  overwrite=TRUE)

file.copy(file.path(repDir,"NtFastasSummary.txt"),
          file.path(resultsDir,"NtFastasSummary.txt"),
		  overwrite=TRUE)
