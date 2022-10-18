
library(ShortRead)
library(Biostrings)
library(stringr)
library(data.table)

fn.fastq <- function(flnm,ln)
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
    ###  Longitudes de secuencia
    sqln <- width(sqq)	
    all.ln <- c(all.ln,sqln)	
    fnm.l <- fnm.l + quantile(sqln,p=c(0.05,0.25,0.5,0.75,0.95),
	                     na.rm=TRUE) * nrds[nchk]
  }
  close(strm)
  return(list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln))
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
flashFiltDir <- "./flashFilt"

###  Lectura de la estructura de las descripciones de las muestras
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Lista de pools en muestras
pools <- unique(samples$Pool.Nm)
###  Lectura de las descripciones de los adaptadores
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Longitud del amplicón mayor en cada pool
max.len.in.pool <- function(p)
{ idx <- which(samples$Pool.Nm==p)
  if(length(idx)==0) return(0)
  idx <- which(primers$Ampl.Nm %in% samples$Primer.ID[idx])
  max(primers$RV.pos[idx]-primers$FW.pos[idx]+1)  
}
pln <- sapply(pools,max.len.in.pool)					  

###  Ficheros resultantes de Flash        
flnms <- list.files(flashFiltDir)
snms <- sub("_flashFilt\\.fastq$","",flnms)

###  Sincronizar
pln <- pln[snms]

###  Loop sobre pools
binsz <- 10
for(i in 1:length(snms))
{ 
  pdf.flnm <- paste("PoolFiltQCbyPos",snms[i],"pdf",sep=".")
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)
  lst <- fn.fastq(file.path(flashFiltDir,flnms[i]),pln[i])
  fvnm <- lst$fvnq
  plot.flash(fvnm)

  nc <- ncol(fvnm)
  fvnm <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm[,iwin:(iwin+9),drop=FALSE]))
  plot.flash(fvnm)

  par(mfrow=c(1,2))  
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
