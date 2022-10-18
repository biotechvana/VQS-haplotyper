
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
  while(length(sqq <- yield(strm)))
  { ###  Longitudes de secuencia
    all.ln <- c(all.ln,width(sqq))	
  }
  close(strm)
  return(all.ln)
}  

#-------------------------------------------------------------------------#

runDir <- "./run"
flashDir <- "./flash"
repDir <- "./reports"
dataDir <- "./data"

###  Lectura de la estructura de las descripciones de las muestras
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Lista de pools en muestras
pools <- unique(samples$Pool.Nm)

###  Ficheros resultantes de Flash        
flnms <- list.files(flashDir)
snms <- sub("_flash\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

pdf(file.path(repDir,"PoolReadLengths.pdf"),paper="a4",width=6,height=10)
par(mfrow=c(2,1),mar=c(5,4,4,2.5)+0.1)
par(mfrow=c(2,1))

###  Loop sobre pools
for(i in 1:length(snms))
{ 
  vln <- fn.fastq(file.path(flashDir,flnms[i]))

  lnfrq <- table(vln)
  x <- as.integer(names(lnfrq))
  y <- as.vector(lnfrq)
  plot(x,y,type="h",xlab="Read length",ylab="Frequency")
  title(main=paste(snms[i],"- read lengths"),line=2)
  par(new=T)
  plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
       axes=F,xlab=NA,ylab=NA,col="blue")
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2,cex.axis=0.8)
  grid()
  idx <- which( y/sum(y) >= 0.05 )
  x[idx]
  title(main=paste("Peaks at",paste(x[idx],collapse=", ")),
        cex.main=0.8,line=0.8)
}

dev.off()
