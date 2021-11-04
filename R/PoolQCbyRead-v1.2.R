
library(ShortRead)
library(Biostrings)
library(stringr)
library(data.table)


fn.fastq <- function(flnm,ln=301)
{
  nrds <- numeric()
  all.ln <- integer()
  all.nl30 <- integer()
  ###  Set streamer on fastq file
  strm <- FastqStreamer(flnm,n=1e6)
  ###  Load fastq file by chuncks
  nchk <- 0
  while(length(sqq <- yield(strm)))
  { nchk <- nchk+1
    nrds[nchk] <- length(sqq)
	###  Phred scores. Codi ASCII-33
    phrsc <- as(quality(sqq),"matrix")
    ###  Bases per sota de Q30    
    nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE))
    all.nl30 <- c(all.nl30,nl30)
    ###  Longituds de sequencia
    sqln <- width(sqq)	
    all.ln <- c(all.ln,sqln)	
  }
  close(strm)
  return(list(all.ln=all.ln,all.nl30=all.nl30))
}  

#-------------------------------------------------------------------------#

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

###  Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
snms <- sub("_flash\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

###  Loop sobre pools
for(i in 1:length(snms))
{ 
  lst1 <- fn.fastq(file.path(flashDir,flnms[i]))

  pdf.flnm <- paste("PoolQCbyRead_",parts[i,"PatID"],".pdf",sep="")
  pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=10.5)
  par(mfrow=c(2,1))
  nl30 <- lst1$all.nl30[lst1$all.nl30<quantile(lst1$all.nl30,p=0.99)]
  hdt <- hist(nl30,breaks=30,main="",xlab="# bases below Q30")
  qs <- quantile(lst1$all.nl30,p=c(0.5,0.8,0.95))
  arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
  arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
  abline(v=qs[1],lty=4,col="blue",lwd=2)
  q95 <- paste("q95% ",qs[3])
  q80 <- paste("q80% ",qs[2])
  med <- paste("q50% ",qs[1])
  text(x=max(hdt$mids)*0.95,y=max(hdt$counts)*0.9,adj=1,cex=0.8,
       paste(q95,q80,med,sep="\n"))
  title(main=parts[i,"PatID"],line=1)
  
  all.fnl30 <- lst1$all.nl30/lst1$all.ln
  fnl30 <- all.fnl30[all.fnl30<quantile(all.fnl30,p=0.99)]
  hdt <- hist(fnl30,breaks=30,main="",
              xlab="fraction of bases below Q30 by read")
  qs <- round(quantile(all.fnl30,p=c(0.5,0.8,0.95)),3)
  arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
  arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
  abline(v=qs[1],lty=4,col="blue",lwd=2)
  q95 <- paste("q95%",qs[3])
  q80 <- paste("q80%",qs[2])
  med <- paste("q50%",qs[1])
  vals <- sprintf("%5.3f",c(qs[3],qs[2],qs[1]))
  vals <- paste(c("q95%","q80%","q50%"),vals,collapse="\n")
  text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.95,adj=1,cex=0.8,vals)
  abline(v=c(0.01,0.02,0.05),lty=4,col="gray")
  p1pct <- sum(all.fnl30<=0.01)/length(all.fnl30)*100
  p2pct <- sum(all.fnl30<=0.02)/length(all.fnl30)*100
  p5pct <- sum(all.fnl30<=0.05)/length(all.fnl30)*100
  vals <- sprintf("%5.1f",c(p1pct,p2pct,p5pct))
  txt <- paste("Pr(f<=",c(1,2,5),"%)  ",vals,"%",sep="")
  txt <- paste(txt,collapse="\n")
  text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.75,adj=1,cex=0.8,txt)
  title(main=parts[i,"PatID"],line=1)
  dev.off()
}
