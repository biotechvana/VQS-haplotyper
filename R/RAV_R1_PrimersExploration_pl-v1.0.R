
library(Biostrings)
library(ShortRead)

repDir <- "./reports"
runDir <- "./run"
dataDir <- "./data"

proj.nm <- "HMar04_RAV"

###  Llegim l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
###  Llegim fitxers globals de HCV
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                    stringsAsFactors=F)
###  Localitzar amplicons de l'experiment
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
samples$Primer.ID <- idx
###  Analitzar consistència en la descripció de mostres
if(sum(is.na(idx)))
  stop("Identificador de primer desconegut")
pnms <- unique(samples$Patient.ID)  

###  Llista de fitxers disponibles
flnms <- list.files(runDir)
flnms <- flnms[ grep("_R1_",flnms) ]
parts <- t(sapply(flnms,function(str) strsplit(str,split="_")[[1]]))
parts <- parts[,-c(3,5)]
colnames(parts) <- c("PatID","SmplID","Read")
snms <- paste(parts[,1],parts[,2],sep="_")
rownames(parts) <- snms

###  Comprobar que tots els fitxers estiguin descrits
if( !all( snms %in% pnms ) )
{ missing.d <- paste(setdiff(snms,pnms),collapse=", ")
  stop(paste("Missing descriptors for fastq files:",missing.d))
}

Ns <- nrow(samples)
Np <- length(snms)
p.cv <- integer(Np)
pr.res <- data.frame(PatID=character(Ns),PrimerID=character(Ns),
                     treads=integer(Ns),
					 fw.match=integer(Ns),rv.match=integer(Ns),
					 stringsAsFactors=FALSE)
k <- 0
for(i in 1:length(snms))
{
  ###  Load R1 fastq file, subset to first 50 nts.
  sqq <- readFastq(dirPath=runDir,patt=flnms[i])
  seqs <- sread(sqq)
  seqs <- seqs[ width(seqs) > 200 ]
  seqs <- subseq( seqs, 1, 50 )
  p.cv[i] <- length(seqs)
  
  idx <- which(samples$Patient.ID == snms[i])
  for(j in idx)
  {
    pr.up <- primers$Primer.FW[ samples$Primer.ID[j] ]
    up.matches <- sum( vcountPattern(pattern=pr.up,subject=seqs,
	                   max.mismatch=3,fixed=FALSE) )
    pr.dn <- primers$Primer.RV[ samples$Primer.ID[j] ]
    dn.matches <- sum( vcountPattern(pattern=pr.dn,subject=seqs,
	                   max.mismatch=3,fixed=FALSE) )
    k <- k+1
    pr.res$PatID[k] =  snms[i]
    pr.res$PrimerID[k] =  primers$Ampl.Nm[ samples$Primer.ID[j] ]
    pr.res$treads[k] =  p.cv[i]
    pr.res$fw.match[k] =  up.matches
    pr.res$rv.match[k] =  dn.matches
  }
}

library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

pdf.flnm <- paste(proj.nm,"R1_PrimersExploration.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

res.mat <- pr.res[,4:5]
ymx <- max(res.mat)*1.2
nms <- paste(pr.res$PatID,pr.res$PrimerID)

bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx))
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))

T.reads <- apply(pr.res[,4:5],2, function(x)
              tapply(x,pr.res$PatID,sum))
yield <- round(rowSums(T.reads)/p.cv*100,1)
barplot(yield,border="navy",col="lavender",las=2,ylim=c(0,100),
        names.arg=rownames(T.reads))
abline(h=0)
title(main="Primer matches global yield (%)")
dev.off()

txt.flnm <- paste(proj.nm,"R1_PrimersExploration.txt",sep="_")
sink(file.path(repDir,txt.flnm))
cat("\nTable of reads identified by primer\n\n")
print(pr.res)
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
cat("\nGlobal yield\n\n")
print(yield)
sink()
