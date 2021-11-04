
library(Biostrings)
library(ShortRead)

#-----------------------------------------------------------------------#

###  Formatar sencers emplenant amb 0 per l'esquerra
######################################################
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}

###  Guarda unes seqüencies en format fasta
#############################################
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm)
}

###  Dóna nom als haplotipus en funció del nombre de mutacions respecte
###   la referència (que pot ser la màster) i de la seva frequencia
###   poblacional, i salva les sequencies en un fitxer fasta.
########################################################################
SaveAllHaplotypes <- function(bseqs,nr,flnm)
{
  code <- "Hpl"
  ##  Determinar diferències respecte la màster
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[1])
  nm <- nmismatch(psa)
  tnm <- table(nm)

  ##  Ordenar per nombre de mutacions
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]

  ##  Numero d'ordre dins de cada nombre de mutacions
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ##  Ordenar per frequencia descendent dins de cada nombre de mutacions
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }

  ##  Calcular frequencia relativa
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotipus
  nms <- paste(code,nm,zeroFillInt2Char(isq,4),sep=".")

  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ##  Salvar a fasta
  write.fasta(DNAStringSet(bseqs),file.path(trimDir,flnm))

  list(bseqs=bseqs,nr=nr,nm=nm)
}

#-----------------------------------------------------------------------#

###  Llegim l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
###  Llegim fitxers globals de HCV
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                    stringsAsFactors=F)

###  Llista de fitxers disponibles
flnms <- list.files(flashDir)
snms <- sub("_flash\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")


###  Inicialitzacions
Ns <- nrow(samples)
Np <- length(flnms)
p.cv <- integer(Np)
p.ok <- integer(Np)
Ns <- Ns*2					 
pr.res <- matrix(0,nrow=Ns,ncol=4)
colnames(pr.res) <-  c("Tot.reads","matches","shorts","fn.reads")
					 
FlTbl <- data.frame(File.Name=character(Ns),Pat.ID=character(Ns),
                    Ampl.Nm=character(Ns),Pr.ID=integer(Ns),
					Str=character(Ns),Pos=integer(Ns),Len=integer(Ns),
					Reads=integer(Ns),Hpls=integer(Ns),
					stringsAsFactors=FALSE)
					
###  Loop sobre pools en samples
k <- 0
pools <- unique(samples$Pool.Nm)
names(p.cv) <- names(p.ok) <- pools
for(i in 1:length(pools))
{
  ###  Identificar fastq del pool
  ip <- which(parts[,"PatID"]==pools[i])
  if(length(ip)==0) next
  ###  Identificar mostres en pool
  idx <- which(samples$Pool.Nm==pools[i])
  if(length(idx)==0) next
  
  ###  Load fastq file
  sqq <- readFastq(dirPath=flashDir,patt=flnms[ip])
  seqs <- sread(sqq)
  seqs <- seqs[width(seqs)>=80] # Eliminar les massa curtes
  p.cv[i] <- length(seqs)
  
  ###  Loop sobre mostres en pool
  for(j in idx)
  {
    ###  primer up matches 
    #################################
    ipr <- grep(samples$Primer.ID[j],primers$Ampl.Nm)
    pr.up <- primers$Primer.FW[ ipr ]
    up.matches <- vmatchPattern(pattern=pr.up,
                                subject=subseq(seqs,start=1,end=80),
                                max.mismatch=max.prdif,fixed=FALSE)
    flags <- elementLengths(up.matches)>=1
    nr <- integer()
	shorts <- integer()
	
    up.flnm <- paste(samples$Patient.ID[j],primers$Ampl.Nm[ipr],"PrFW.fna",
	                 sep=".")
    if(sum(flags))
    { ###  Aquest pas és delicat, startIndex() té NULLs a la llista.
      pos <- sapply(startIndex(up.matches)[flags],function(x) x[[1]])
      ###  Trim up seqs
      seqs.up <- seqs[flags]
      st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
	  trim.len <- primers$RV.tpos[ipr]-primers$FW.tpos[ipr]+1
      last <- pmin(st+trim.len-1,width(seqs.up))
      seqs.up <- subseq(seqs.up,start=st,end=last)
      ###  Exclude sequences not covering the full amplicon
      shorts <- width(seqs.up)<trim.len
      seqs.up <- seqs.up[!shorts]

      if(length(seqs.up))
      { ###  Collapse sequences to haplotypes+frequencies
        sqtbl <- sort(table(as.character(seqs.up)),decreasing=TRUE)
        bseqs <- names(sqtbl)                
        names(bseqs) <- 1:length(bseqs)
        nr <- as.integer(sqtbl)
        lst <- SaveAllHaplotypes(bseqs,nr,up.flnm)
		nr <- lst$nr
      }
    }
	k <- k+1
	FlTbl$File.Name[k] <- up.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[j]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "fw"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- trim.len
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   

	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  
    ###  primer dn matches 
    #################################
    pr.dn <- primers$Primer.RV[ ipr ]
    dn.matches <- vmatchPattern(pattern=pr.dn,
                                subject=subseq(seqs,start=1,end=80),
                                max.mismatch=max.prdif,fixed=FALSE)
    flags <- elementLengths(dn.matches)>=1
    nr <- integer()
    shorts <- integer()
	
    dn.flnm <- paste(samples$Patient.ID[j],primers$Ampl.Nm[ipr],"PrRV.fna",
	                 sep=".")
    if(sum(flags))
    { ###  Aquest pas és delicat, startIndex() té NULLs a la llista.
      pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])

      ###  Reverse complement down matches
      seqs.dn <- reverseComplement(seqs[flags])
      ###  ... and match up primer
      dn.matches <- vmatchPattern(pattern=pr.up,
                                  subject=subseq(seqs.dn,start=1,end=80),
                                  max.mismatch=max.prdif,fixed=FALSE)
      flags <- elementLengths(dn.matches)>=1
      if(sum(flags))
	  { ###  Trim reversed dn seqs
        seqs.dn <- seqs.dn[flags]
        pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])
        st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
	    trim.len <- primers$RV.tpos[ipr]-primers$FW.tpos[ipr]+1
        last <- pmin(st+trim.len-1,width(seqs.dn))
        seqs.dn <- subseq(seqs.dn,start=st,end=last)
        ###  Exclude sequences not covering the full amplicon
        shorts <- width(seqs.dn)<trim.len
        seqs.dn <- seqs.dn[!shorts]

        if(length(seqs.dn))
        { ###  Collapse sequences to haplotypes+frequencies
          sqtbl <- sort(table(as.character(seqs.dn)),decreasing=TRUE)
          bseqs <- names(sqtbl)                
          names(bseqs) <- 1:length(bseqs)
          nr <- as.integer(sqtbl)
	      lst <- SaveAllHaplotypes(bseqs,nr,dn.flnm)
          nr <- lst$nr
		}
	  }
	}  

    k <- k+1
	FlTbl$File.Name[k] <- dn.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[j]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "rv"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- trim.len
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   
	
	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  }
}


###  Plot results
library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

fw.idx <- which(FlTbl$Str=="fw")
rv.idx <- which(FlTbl$Str=="rv")
mprres <- data.frame(PatID=FlTbl$Pat.ID[fw.idx],
                     PrimerID=FlTbl$Ampl.Nm[fw.idx],
                     Treads=pr.res[fw.idx,1],
					 Shorts=pr.res[fw.idx,3]+pr.res[rv.idx,3],
                     FW.match=pr.res[fw.idx,4],
                     RV.match=pr.res[rv.idx,4],
					 Fn.reads=pr.res[fw.idx,4]+pr.res[rv.idx,4],
                     stringsAsFactors=FALSE)
mres <- mprres
T.reads <- apply(mres[,5:6],2, function(x)
              tapply(x,mres$PatID,sum))

pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

ymx <- max(rowSums(T.reads))*1.2
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by sample (# reads)")

res.mat <- mres[,5:6]
ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))

yield <- round(p.ok/p.cv*100,1)
barplot(yield,border="navy",col="lavender",las=2,ylim=c(0,100),
        names.arg=names(p.cv))
abline(h=0)
title(main="Primer matches by pool (%)")
dev.off()

PoolTbl <- data.frame(FlashReads=p.cv,PrimerReads=p.ok,Pct=yield)

txt.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.txt",sep="_")
sink(file.path(repDir,txt.flnm))
cat("\nTable of reads identified by primer\n\n")
print(mprres)
cat("\n")
print(FlTbl[,-1])
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
cat("\nYield by pool\n\n")
print(PoolTbl)
sink()

###  Save tables

# FlTbl <- FlTbl[FlTbl$Reads>0, ]  # Erase rows of null files
save(FlTbl,PoolTbl,pr.res,file=file.path(repDir,"SplittedReadsFileTable.RData"))
sink(file.path(repDir,"SplittedReadsFileTable.txt"))
print(FlTbl)
sink()


###  Plot on A4 horizontally
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash-hz.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+0.1)

ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
pgrp <- factor(mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
points(colMeans(bp),rep(min(res.mat)/2,nrow(res.mat)),
       pch=(1:nlevels(pgrp))[pgrp],cex=0.8,font=2)
	   
par(mfrow=c(1,2))
ymx <- max(rowSums(T.reads))*1.2
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

barplot(yield,border="navy",col="lavender",las=2,ylim=c(0,100),
        names.arg=rownames(p.cv))
abline(h=0)
title(main="Primer matches by pool (%)")

ymx <- max(c(res.mat[,1],res.mat[,2]))
boxplot(res.mat[,1]~samples$Region,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(samples$Region)),a=0.15),res.mat[,1],
       pch="+",cex=0.8)
title(main="Primers identified on forward reads")	   

boxplot(res.mat[,2]~samples$Region,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(samples$Region)),a=0.15),res.mat[,2],
       pch="+",cex=0.8)
title(main="Primers identified on reverse reads")	   

dev.off()
