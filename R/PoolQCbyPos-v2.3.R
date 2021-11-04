
library(ShortRead)
library(Biostrings)
library(stringr)
library(data.table)

fn.fastq <- function(flnm,ln=301)
{
  fnm.q <- matrix(0,nrow=5,ncol=ln)
  fnm.l <- numeric(5)
  nrds <- numeric()
  all.ln <- integer()
  
  ###  Set streamer on fastq file
  strm <- FastqStreamer(flnm,n=1e6)
  ###  Load fastq file by chuncks
  nchk <- 0
  while(length(sqq <- yield(strm)))
  { nchk <- nchk+1
    nrds[nchk] <- length(sqq)
	###  Phred scores. Codi ASCII-33
    phrsc <- as(quality(sqq),"matrix")
	nc <- min(ln,ncol(phrsc))
    fnm.q[,1:nc] <- fnm.q[,1:nc] + 
	               apply(phrsc,2,quantile,p=c(0.05,0.25,0.5,0.75,0.95),
	                     na.rm=TRUE)[,1:nc] * nrds[nchk]
    ###  Longituds de sequencia
    sqln <- width(sqq)	
    all.ln <- c(all.ln,sqln)	
    fnm.l <- fnm.l + quantile(sqln,p=c(0.05,0.25,0.5,0.75,0.95),
	                     na.rm=TRUE) * nrds[nchk]
  }
  close(strm)
  return(list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln))
}  

plot.R1R2 <- function(fvnm1,fvnm2,snm,SW=FALSE)
{ nc <- ncol(fvnm1)
  plot(1:nc,fvnm1[2,],type="l",col="navy",ylim=c(0,40),xaxt="n",
       xlab="",ylab="Phred score")
  lines(1:nc,fvnm1[3,],type="l",lwd=2,col="navy")
  lines(1:nc,fvnm1[4,],type="l",col="navy")
  lines(1:nc,fvnm1[1,],type="l",lty=3,col="navy")
  lines(1:nc,fvnm1[5,],type="l",lty=3,col="navy")
  axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
  nc <- ncol(fvnm2)
  lines(1:nc,fvnm2[2,],type="l",col="maroon")
  lines(1:nc,fvnm2[3,],type="l",lwd=2,col="maroon")
  lines(1:nc,fvnm2[4,],type="l",col="maroon")
  lines(1:nc,fvnm2[1,],type="l",lty=3,col="maroon")
  lines(1:nc,fvnm2[5,],type="l",lty=3,col="maroon")
  abline(h=c(20,23,30,33),lty=4,col="gray")
  legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
         legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
  legend("bottom",lwd=2,col=c("navy","maroon"),legend=c("R1","R2"),horiz=TRUE)		 
  if(SW)
    title(main=paste("Quality profile by SW (size 10,step 1):",snm))
  else
    title(main=paste("Quality profile by position:",snm))
}

plot.flash <- function(fvnm)
{ nc <- ncol(fvnm)
  bxp.dt <- list(stats=fvnm,names=1:ncol(fvnm))
  plot(1:nc,fvnm[2,],type="l",col="darkgreen",ylim=c(0,40),xaxt="n",
       xlab="Position",ylab="Phred score")
  lines(1:nc,fvnm[3,],type="l",lwd=2,col="darkgreen")
  lines(1:nc,fvnm[4,],type="l",col="darkgreen")
  lines(1:nc,fvnm[1,],type="l",lty=3,col="darkgreen")
  lines(1:nc,fvnm[5,],type="l",lty=3,col="darkgreen")
  axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
  abline(h=c(20,23,30,33),lty=4,col="gray")
  legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
         legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
  legend("bottom",lwd=2,col="darkgreen",legend="Flash reads",horiz=TRUE)		 
}	

#-------------------------------------------------------------------------#

runDir <- "./run"
flashDir <- "./flash"
repDir <- "./reports"
dataDir <- "./data"

###  Llegim l'estructura de descripció de mostres
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Llista de pools en samples
pools <- unique(samples$Pool.Nm)
###  Llegim descriptors de primers
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Longitud de l'amplicó major en cada pool
max.len.in.pool <- function(p)
{ idx <- which(samples$Pool.Nm==p)
  if(length(idx)==0) return(0)
  idx <- which(primers$Ampl.Nm %in% samples$Primer.ID[idx])
  max(primers$RV.pos[idx]-primers$FW.pos[idx]+1)  
}
pln <- sapply(pools,max.len.in.pool)					  
###  Amb MID+M13
pln <- pln + 2*(20+10)

###  Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
snms <- sub("_flash\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

###  Fitxers R1 i R2 origen
R1.flnms <- paste(snms,"_L001_R1_001.fastq.gz",sep="")
R2.flnms <- paste(snms,"_L001_R2_001.fastq.gz",sep="")

###  Sincronitzar
pln <- pln[parts[,"PatID"]]

###  Loop sobre pools
binsz <- 10
for(i in 1:length(snms))
{ 
  pdf.flnm <- paste("PoolQCbyPos",parts[i,1],"pdf",sep=".")
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)
  lst1 <- fn.fastq(file.path(runDir,R1.flnms[i]))
  fvnm1 <- lst1$fvnq
  lst2 <- fn.fastq(file.path(runDir,R2.flnms[i]))
  fvnm2 <- lst2$fvnq
  plot.R1R2(fvnm1,fvnm2,snms[i])
  lst <- fn.fastq(file.path(flashDir,flnms[i]),pln[i])
  fvnm <- lst$fvnq
  plot.flash(fvnm)

  nc <- ncol(fvnm1)
  fvnm1 <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm1[,iwin:(iwin+9),drop=FALSE]))
  nc <- ncol(fvnm2)
  fvnm2 <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm2[,iwin:(iwin+9),drop=FALSE]))
  plot.R1R2(fvnm1,fvnm2,snms[i],TRUE)
  nc <- ncol(fvnm)
  fvnm <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm[,iwin:(iwin+9),drop=FALSE]))
  plot.flash(fvnm)

  par(mfrow=c(1,2))  
  stats=cbind(lst1$fvnl,lst2$fvnl,lst$fvnl)
  colnames(stats) <- c("R1","R2","Flash")
  bxp.dt <- list(stats=stats,names=colnames(stats))
  bxp(bxp.dt,pars=list(boxfilll="lavender"),border="navy",
      ylab="Read length",las=2,ylim=c(0,max(stats)))
  title(main="Read length distributions")
  
  lnfrq <- table(lst$all.ln)
  x <- as.integer(names(lnfrq))
  y <- as.vector(lnfrq)
  plot(x,y,type="h",xlab="Read length",ylab="Frequency")
  title(main="Read lengths (Flash reads)")
  par(new=T)
  plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
       axes=F,xlab=NA,ylab=NA,col="blue")
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
  grid()
  
  dev.off()
}
