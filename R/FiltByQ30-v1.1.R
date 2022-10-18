
library(ShortRead)
library(Biostrings)
library(stringr)
library(data.table)

###  Lectura de la estructura de las descripciones de las muestras
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Lista de pools en muestras
pools <- unique(samples$Pool.Nm)
###  Lectura de las descripciones de los adaptadores
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)

###  Ficheros resultantes de Flash        
flnms <- list.files(flashDir)
snms <- sub("_flash\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

oflnms <- paste(parts[,"PatID"],"flashFilt.fastq",sep="_")
oflnms <- file.path(flashFiltDir,oflnms)
flnms <- file.path(flashDir,flnms)

###  Cargar el rendimiento de flash
load(file.path(repDir,"FLASH_table.RData"))

###  Loop sobre pools
freads <- integer(length(snms))
for(i in 1:length(snms))
{ 
  raw.rds <- filt.rds <- 0
  
  ###  Set streamer on fastq file
  strm <- FastqStreamer(flnms[i],n=1e6)
  ###  Load fastq file by chuncks
  if(file.exists(oflnms[i]))
    file.remove(oflnms[i])
  appfl <- "w"
  while(length(sqq <- yield(strm)))
  { raw.rds <- raw.rds+length(sqq)
    ###  Phred scores. Codi ASCII-33
    phrsc <- as(quality(sqq),"matrix")
    ###  Bases por debajo de Q30    
    nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE))
    ###  Longitudes de secuencia
    sqln <- width(sqq)	
    ###  Fracciones
	fnl30 <- nl30/sqln
	###  Filtro
	sqq <- sqq[fnl30<=ThrQ30]
	filt.rds <- filt.rds+length(sqq)
	writeFastq(sqq,oflnms[i],mode=appfl,compress=TRUE)
	appfl <- "a"
  }
  close(strm)
  freads[i] <- filt.rds  
}
flash.res$FiltQ30 <- freads
flash.res$Raw <- flash.res$Extended+flash.res$NoExtd
save(flash.res,file=file.path(repDir,"FLASH_table.RData"))

res <- data.matrix(flash.res[,c("Raw","Extended","FiltQ30")])

pdf(file.path(repDir,"FiltQ30.barplot.pdf"),paper="a4",width=6.6,height=10)
par(mfrow=c(2,1))
library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
ymx <- max(res)*1.15
bp <- barplot(t(res),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
               xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))

resy <- round(res/res[,1]*100,1)
ymx <- 115
bp <- barplot(t(resy),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
               xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))
y <- min(resy)/2
text(x=as.vector(bp),y=y,lab=as.vector(t(resy)),cex=0.5,font=2,srt=90,
     col=pal1[1:3])
dev.off()

sink(file.path(repDir,"FiltQ30_report.txt"))
print(flash.res)
sink()

