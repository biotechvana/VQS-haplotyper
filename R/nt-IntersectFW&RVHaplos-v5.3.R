####################################################
###        HAPLOTIPOS de CONSENSO FW y RV        ###
####################################################
options(max.print=999999)
source(file.path(codeDir,"global.v4.6.R"))
library(data.table)

###  Formatear enteros rellenando con 0 por la izquierda
######################################################
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}

###  Guarda secuencias en formato fasta
#############################################
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm)
}

###  Alinear las distribuciones de haplotipos de dos poblaciones
################################################################
PopsAlgnHist <- function(IDsA,seqsA,IDsB,seqsB)
{ nA <- IDsA$nseqs
  names(nA) <- seqsA
  nB <- IDsB$nseqs
  names(nB) <- seqsB
  nms <- union(seqsA,seqsB)    # todos: A + B, elimina las secuencias duplicadas
  nb <- length(nms)
  pA <- integer(nb)
  idx <- which(nms %in% seqsA) # los de A en todos
  pA[idx] <- nA[nms[idx]]  
  pB <- integer(nb)
  idx <- which(nms %in% seqsB) # los de B en todos
  pB[idx] <- nB[nms[idx]]  
  list(pA=pA,pB=pB,Hpl=nms)
}

###  Confrontar haplotipos en FW y RV
Intersect.FWRV <- function(IDsA,seqsA,IDsB,seqsB)  
{
  ###  Alinear distribuciones
  lst <- PopsAlgnHist(IDsA,seqsA,IDsB,seqsB)
  ###  Haplotipos comunes y solapamiento global
  fl <- lst$pA>0 & lst$pB>0
  ov.a <- sum(lst$pA[fl]+lst$pB[fl])/sum(lst$pA+lst$pB)
  ###  Frecuencias y solapamiento por intersección
  pFW <- (lst$pA/sum(lst$pA))[fl]
  pRV <- (lst$pB/sum(lst$pB))[fl]
  p <- pmin(pFW,pRV) # intersección
  ov.i <- sum(p)  
  p <- p/sum(p)      # renormalización

  list(p=p,seqs=lst$Hpl,ov.i=ov.i,ov.a=ov.a,pA=lst$pA,pB=lst$pB)
}

### Representa las distribuciones alineadas de los haplotipos
##############################################################
PlotHplHistos <- function(tt,pA,pB,p)
{ pA <- pA/sum(pA)
  pB <- pB/sum(pB)
  ymx <- max(c(pA,pB))
  barplot(pA,ylim=c(0,ymx)); abline(h=0)
  title(main="FW strand haplotypes histogram",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  barplot(pB,ylim=c(0,ymx)); abline(h=0)
  title(main="RV strand haplotypes histogram",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  barplot(p); abline(h=0)
  title(main="Intersected haplotypes histogram",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
}


###  Da nombre a los haplotipos en función del número de mutaciones
###  respecto la referencia (que puede ser la máster) y de su frecuencia
###  poblacional, y guarda las secuencias en un fichero fasta.
########################################################################
SaveHaplotypes <- function(flnm,bseqs,nr)
{
  ##  Determinar diferencias respecto a la máster
  i.mast <- which.max(nr)
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[i.mast])
  nm <- nmismatch(psa)
  tnm <- table(nm)

  ##  Ordenar por número de mutaciones
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]

  ##  Número de orden dentro de cada número de mutaciones
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ##  Ordenar por frecuencias descendente dentro de cada número de mutaciones
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }

  ##  Calcular frecuencia relativa
  frq <- round(nr/sum(nr)*100,2)

  ## Nombre completo por cada haplotipo
  nms <- paste("Hpl",nm,zeroFillInt2Char(isq,4),sep=".")

  ##  Cabecera fasta con nombre, núemro de reads y frecuencia relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ##  Guardar a fasta
  write.fasta(DNAStringSet(bseqs),flnm)

  list(bseqs=bseqs,nr=nr,nm=nm)
}

###  Filtrar por encima de un nivel de abundancia (0 por defecto para reporte todos los haplotipos)
filter.haplos <- function(lst,pct.cut,min.rd)
{ nr <- lst$IDs$nseqs
  pct <- nr/sum(nr)*100
  fl.in <- (pct >= pct.cut & nr >= min.rd)
  lst$IDs <- lst$IDs[fl.in,]
  lst$seqs <- lst$seqs[fl.in]
  lst
}

SaveNotCommon <- function(lst,desc,thr)
{ lost <- character()
  ###  dominante
  fl <- lst$pA>0 & lst$pB>0
  sq <- lst$seqs[fl]
  nr <- lst$pA[fl]+lst$pB[fl]
  idx <- which.max(nr)
  sq <- sq[idx]
  p <- round(nr[idx]/sum(nr)*100,2)
  names(sq) <- paste("Master.",nr[idx],"|",p,sep="")
  ###  forward unicos bastante abundantes
  p  <- lst$pA/sum(lst$pA)*100
  fl <- lst$pA>0 & lst$pB==0 & p>=thr
  if( sum(fl) )
  { lseqs <- lst$seqs[fl]
    names(lseqs) <- paste("Lost.fw.",lst$pA[fl],"|",round(p[fl],2),sep="")
    lost <- c(lost,lseqs)
  }
  ###  reverse unicos bastante abundantes
  p  <- lst$pB/sum(lst$pB)*100
  fl <- lst$pA==0 & lst$pB>0 & p>=thr
  if( sum(fl) )
  { lseqs <- lst$seqs[fl]
    names(lseqs) <- paste("Lost.rv.",lst$pB[fl],"|",round(p[fl],2),sep="")
    lost <- c(lost,lseqs)
  }

  ##  Guardar en fasta
  if(length(lost))
  { lost <- c(sq,lost)
    flnm <- paste(desc$Pat.ID[1],desc$Ampl.Nm,"lost.fna",sep=".")
    flnm <- file.path(joinDir,flnm)
    write.fasta(DNAStringSet(lost),flnm)
  }
}


###   MAIN LINE
#############################################################################
library(Biostrings)
library(stringr)         
options(width=80)
library(RColorBrewer)
cls <- brewer.pal(8,"Dark2")

source(file.path(codeDir,"global.v4.6.R"))

###  Lectura de la estructura de las descripciones de las muestras
###################################################################
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                    stringsAsFactors=F)
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                    stringsAsFactors=F)
RefSeqs <- as.character(readDNAStringSet(file.path(dataDir,"AmpliconRefSeqs.fna")))
idx <- sapply(samples$RefSeq.ID,function(pr) which(names(RefSeqs)==pr)[1])
samples$Primer.idx <- idx

##  Leer la tabla de ficheros a tratar y generar nombres
#########################################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))

###  Eliminar muestras sin pareja FW o RV
if(sum(FlTbl$Str=="fw") != sum(FlTbl$Str=="rv"))
{ idx.fw <- which(FlTbl$Str=="fw")
  idx.rv <- which(FlTbl$Str=="rv")
  noc <- c(setdiff(idx.fw,idx.rv-1),setdiff(idx.rv-1,idx.fw))
  FlTbl <- FlTbl[-noc,]
}

flnms <- FlTbl$File.Name
nms.fw <- flnms[ FlTbl$Str=="fw" ]
nms.rv <- flnms[ FlTbl$Str=="rv" ]

flnms <- data.frame(fw=nms.fw, rv=nms.rv,
                    stringsAsFactors=FALSE )
flnms$fw <- file.path(filtDir,flnms$fw)
flnms$rv <- file.path(filtDir,flnms$rv)

out.flnms <- sub("\\.PrFW\\.","\\.",nms.fw)
out.flnms <- file.path(joinDir,out.flnms)


###  Inicio del informe
########################
sink(file=file.path(repDir,"FW+RV.10.Intersects.txt"))

cat("\n   FW + RV HAPLOTYPES INTERSECTIONS")
cat("\n======================================\n")
cat("\nCutting FW and RV at ",a.cut,"% followed by haplotypes intersection.",
    sep="")
if(method=="Sum")
  cat("\nFrequencies as sum of reads FW+RV.\n\n")
if(method=="Intersect")
  cat("\nFrequencies estimated from intersection and renormalization.\n\n")

pdf.flnm <- "FW+RV.10.HaplosBarPlots.pdf"
pdf(file=file.path(repDir,pdf.flnm),paper="a4",width=7,height=11) 
par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
par(mfrow=c(3,1))


rlst <- list(nrow(flnms))
n <- nrow(flnms)
rdf.fw <- data.frame(fw.all=integer(n),fw.lowf=integer(n),
                     fw.in=integer(n),fw.unq=integer(n),fw.com=integer(n))
rdf.rv <- data.frame(rv.all=integer(n),rv.lowf=integer(n),
                     rv.in=integer(n),rv.unq=integer(n),rv.com=integer(n))
rdf.gbl <- data.frame(all=integer(n),lowf=integer(n),unq=integer(n),
                      ovrlp=numeric(n),common=numeric(n),Fn.rd=integer(n))
pct.m <- numeric(n)
FlTbl <- FlTbl[FlTbl$Str=="fw",]

for(i in 1:n)
{  pr.idx <- FlTbl$Pr.ID[i]

   if( ! file.exists(flnms$fw[i]) ) next
   if( ! file.exists(flnms$rv[i]) ) next

   cat("\nProcessing files:")
   cat("\n   ",flnms$fw[i])
   cat("\n   ",flnms$rv[i],"\n")
   cat("\n      Positions ",primers$FW.tpos[pr.idx]," to ",
       primers$RV.tpos[pr.idx],"\n",sep="")
   if(!is.null(primers$FW.tpos))
   { off <- primers$FW.tpos[pr.idx]-1
   } else {
     off <- primers$FW.pos[pr.idx]+nchar(primers$Primer.FW[pr.idx])
   }

   ###  Read and intersect
   fw.lst <- read.ampl.seqs(flnms$fw[i],mnr=1)
   anrd <- rdf.fw$fw.all[i] <- sum(fw.lst$IDs$nseqs)
   fw.lst <- filter.haplos(fw.lst,a.cut,min.rd)
   rdf.fw$fw.lowf[i] <- anrd-sum(fw.lst$IDs$nseqs)

   rv.lst <- read.ampl.seqs(flnms$rv[i],mnr=1)
   anrd <- rdf.rv$rv.all[i] <- sum(rv.lst$IDs$nseqs)
   rv.lst <- filter.haplos(rv.lst,a.cut,min.rd)
   rdf.rv$rv.lowf[i] <- anrd-sum(rv.lst$IDs$nseqs)

   lst <- Intersect.FWRV( fw.lst$IDs,fw.lst$seqs,rv.lst$IDs,rv.lst$seqs)

   cat("\n  Haplotypes distribution overlap: ",round(lst$ov.i*100,2),"%",
       sep="")
   cat("\n       Overlap as reads in common: ",round(lst$ov.a*100,2),"%\n",
       sep="")
   cat("\n       Input reads: ",sum(lst$pA),", ",sum(lst$pB),sep="")
   fl <- lst$pA>0 & lst$pB>0  ###  Haplos comuns
   cat("\n   Reads in common: ",sum(lst$pA[fl]),", ",sum(lst$pB[fl]),sep="")
   cat("\n        Reads lost: ",sum(lst$pA[!fl]),", ",sum(lst$pB[!fl]),sep="")
   cat("\n\n")

   rdf.fw$fw.in[i]  <- sum(lst$pA)
   rdf.rv$rv.in[i]  <- sum(lst$pB)
   rdf.fw$fw.com[i] <- sum(lst$pA[fl])
   rdf.rv$rv.com[i] <- sum(lst$pB[fl])
   rdf.fw$fw.unq[i] <- sum(lst$pA[!fl])
   rdf.rv$rv.unq[i] <- sum(lst$pB[!fl])

   rdf.gbl$all[i]    <- rdf.fw$fw.all[i]+rdf.rv$rv.all[i]
   rdf.gbl$ovrlp[i]  <- round(lst$ov.i*100,2)
   rdf.gbl$common[i] <- round(lst$ov.a*100,2)
   rdf.gbl$Fn.rd[i]  <- rdf.fw$fw.com[i]+rdf.rv$rv.com[i]
   rdf.gbl$lowf[i]   <- rdf.fw$fw.lowf[i] + rdf.rv$rv.lowf[i]
   rdf.gbl$unq[i]    <- rdf.fw$fw.unq[i] + rdf.rv$rv.unq[i]

   ###  No overlap, skip.
   if(sum(fl)==0) 
   { cat("\n--------------------------------------------------\n")
     next
   }
               
   ###  Plot aligned haplotypes bar plots
   tt <- paste(FlTbl$Pat.ID[i],FlTbl$PAmpl.Nm[i],sep=".")
   PlotHplHistos(tt,lst$pA,lst$pB,lst$p)

   ###  Simple suma de reads en FW y RV, proporción como un promedio ponderado
   if(method=="Sum")        
        rds <- lst$pA[fl]+lst$pB[fl]
   ###  Interseccion+renormalización y estimacion de reads 
   if(method=="Intersect")  
        rds <- round(sum(lst$pA[fl]+lst$pB[fl])*lst$p)

   ###  Save common haplos with estimated frequencies
   lst0 <- lst
   lst <- SaveHaplotypes(out.flnms[i],lst$seqs[fl],rds)

   ###  Save haplos not in FW or in RV, but at abundance above
   SaveNotCommon(lst0,FlTbl[i,],ni.thr)
	 
   cat("Common haplotypes:\n\n")
   print(data.frame(Haplotypes=names(lst$bseqs)))
   cat("\nNumber of haplotypes by mutation number:")
   print(table(lst$nm))
   cat("\nReads by number of mutations:\n")
   print(tapply(lst$nr,lst$nm,sum))

   ###  Report point mutations
   if(length(lst$bseqs)>1)
   { refseq.idx <- which(grepl(FlTbl$Ampl.Nm[i], names(RefSeqs)))
     muts <- SummaryMuts.w(lst$bseqs,lst$nr,off,RefSeqs[refseq.idx][[1]])
     cat("\nObserved point mutations:\n")
     print(muts)
     cat("\n")
   }
   pct.m[i] <- max(round(lst$nr/sum(lst$nr)*100,2))
   cat("\n--------------------------------------------------\n")
}
dev.off()

cat("\n\n   SUMMARY RESULTS BY READS NUMBER")
cat("\n=====================================\n\n")
fl.fw <- FlTbl$Str=="fw"
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.fw,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.rv,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.gbl,stringsAsFactors=FALSE)
print(frdf)

cat("\n\nTotal counts:\n\n")
tots.fw <- apply(rdf.fw,2,sum)
print(tots.fw)
cat("\n")
tots.rv <- apply(rdf.rv,2,sum)
print(tots.rv)
cat("\n")
tots.gbl <- apply(rdf.gbl[,c(1:3,6)],2,sum)
print(tots.gbl)
cat("\n")

cat("\nPercentage:\n\n")
ptts <- round(tots.fw/tots.fw[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.rv/tots.rv[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.gbl/tots.gbl[1]*100,2)
print(ptts)
cat("\n")

sink()

txt.flnm <- "FW+RV.10.Intersects-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm),split=TRUE)

cat("\n   FW + RV HAPLOTYPES INTERSECTIONS")
cat("\n======================================\n")
cat("\nCutting FW and RV at ",a.cut,"% followed by haplotypes intersection.\n",
    sep="")
if(method=="Sum")
  cat("\nFrequencies as sum of reads FW+RV.\n\n")
if(method=="Intersect")
  cat("\nFrequencies estimated from intersection and renormalization.\n\n")

cat("\nSUMMARY RESULTS BY READS NUMBER")
cat("\n===============================\n\n")
fl.fw <- FlTbl$Str=="fw"
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.fw,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.rv,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.gbl,stringsAsFactors=FALSE)
print(frdf)

cat("\n\nTotal counts:\n\n")
tots.fw <- apply(rdf.fw,2,sum)
print(tots.fw)
cat("\n")
tots.rv <- apply(rdf.rv,2,sum)
print(tots.rv)
cat("\n")
tots.gbl <- apply(rdf.gbl[,c(1:3,6)],2,sum)
print(tots.gbl)
cat("\n")

cat("\nPercentage:\n\n")
ptts <- round(tots.fw/tots.fw[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.rv/tots.rv[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.gbl/tots.gbl[1]*100,2)
print(ptts)
cat("\n")

sink()


pdf.flnm2 <- "IntersectBarplots.pdf"
pdf(file.path(repDir,pdf.flnm2),paper="a4",width=5.5,height=10)
par(mfrow=c(2,1))
pal <- cls[c(1,3,2)]
mbp <- data.matrix(frdf[,c(8,5,4)])
mbp.nms <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
omar <- par(mar=c(6,3.5,3,2))
barplot(t(mbp),col=pal,ylim=c(0,max(rowSums(mbp))*1.2),names.arg=mbp.nms,
    las=2,cex.names=0.7)
legend("top",horiz=TRUE,fill=pal,legend=c("Common","Unique","Low freq."),
    cex=0.8)
par(mar=omar)
barplot(frdf$common,col="lavender",ylim=c(0,100),las=2)
title(main="FW + RV intersection yield",line=1)

all.nms <- mbp.nms
ytbl <- cbind(All=frdf$all,Passed=frdf$Fn.rd)
rownames(ytbl) <- all.nms
barplot(t(ytbl),beside=TRUE,col=cls[1:2],ylim=c(0,max(ytbl)*1.2),las=2,
     names.arg=rownames(ytbl),cex.names=0.7)
legend("top",horiz=TRUE,fill=cls,legend=colnames(ytbl),cex=0.8)

par(mar=omar)
barplot(as.vector(ytbl[,2]/ytbl[,1]*100),col="lavender",ylim=c(0,100),las=2)
title(main="Ab filter & intersecton yield",line=1)
dev.off()
	

###  Global yield by step

load(file=file.path(repDir,"FLASH_table.RData"))
rownames(flash.res) <- str_extract(rownames(flash.res),"^[A-Za-z0-9\\.-]+")

load(file=file.path(repDir,"SplittedReadsFileTable.RData"))

p.nms <- rownames(flash.res)
rds.raw <- flash.res[,1]+flash.res[,2]
rds.flash <- flash.res[,1]

fltbl.pool <- sapply(1:nrow(FlTbl), function(i)
               samples$Pool.Nm[which(samples$Patient.ID==FlTbl$Pat.ID[i] &
		                             samples$Primer.ID==FlTbl$Ampl.Nm[i])[1] ])
rds.dmult <- tapply(FlTbl$Reads,fltbl.pool,sum)

frdf.pool <- sapply(1:nrow(frdf), function(i)
               samples$Pool.Nm[ which(samples$Patient.ID==frdf$Pat.ID[i] &
		                              samples$Primer.ID==frdf$Ampl.Nm[i])[1] ])
rds.filt <- tapply(frdf$all,frdf.pool,sum)[p.nms]
rds.ints <- tapply(frdf$Fn.rd,frdf.pool,sum)[p.nms]

if(is.null(flash.res$FiltQ30))
{  gbly <- cbind(raw=rds.raw,flash=rds.flash, dmult=rds.dmult,
                 filt=rds.filt,ints=rds.ints)
} else {
   gbly <- cbind(raw=rds.raw,flash=rds.flash, FltQ30=flash.res$FiltQ30,
                 dmult=rds.dmult,filt=rds.filt,ints=rds.ints)
}				 


pdf.flnm3 <- "GlobalYieldBarplots.pdf"
pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
par(mfrow=c(2,1))

pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
m <- ncol(gbly)
ymx <- max(gbly)*1.2
bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
grid(nx=NA,ny=NULL)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

par(mar=c(5,8,4,6))
boxplot(frdf$Fn.rd,border="gray",ylab="# of reads")
points(jitter(rep(1,nrow(frdf)),a=0.10),frdf$Fn.rd,pch="+",cex=0.8)
title(main="Final coverage")

dev.off()
			  
txt.flnm2 <- "GlobalYield-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm2),split=TRUE)

cat("\n   Global yield by analysis step")
cat("\n===================================\n")
cat("\nIn number of reads:\n\n")
gbly <- rbind(gbly,TOTAL=colSums(gbly))
print(gbly)	
		
cat("\nIn percentage referred to raw reads:\n\n")
print(round(gbly/gbly[,1]*100,2))			
sink()

file.copy( file.path(repDir,txt.flnm), file.path(exportDir,txt.flnm),
           overwrite=TRUE )
file.copy( file.path(repDir,txt.flnm2), file.path(exportDir,txt.flnm2),
           overwrite=TRUE )
file.copy( file.path(repDir,pdf.flnm), file.path(exportDir,pdf.flnm),
           overwrite=TRUE )
file.copy( file.path(repDir,pdf.flnm2), file.path(exportDir,pdf.flnm2),
           overwrite=TRUE )
file.copy( file.path(repDir,pdf.flnm3), file.path(exportDir,pdf.flnm3),
           overwrite=TRUE )
