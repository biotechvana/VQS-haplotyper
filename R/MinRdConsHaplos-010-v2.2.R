########################################################
###     FILTER CONSENSUS HAPLOTYPES BY ABUNDANCE     ###
########################################################

options(max.print=999999)

library(Biostrings)
source(file.path(codeDir,"global.v4.6.R"))
library(data.table)


###  Read amplicon aligned sequences 
######################################  
read.ampl.seqs <- function(flnm,mnr=2)
{
  seqs <- as.character(readDNAStringSet(flnm))
  IDstr <- names(seqs)
  n <- length(IDstr)
  nms <- character(length=n)
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1")
  for(j in 1:n)
  { strs <- strsplit(IDstr[j],split="\\|")[[1]]
    nms[j] <- strs[1]
    sts[j,] <- as.numeric(strs[2:3])
  }
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  nall <- nrow(IDs)
  tnr <- sum(IDs$nseqs)
  ###  Filter by minimum reads by haplotype
  flags <- IDs$nseqs >= mnr
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}

###  Select haplos by a given reads threshold
###  Save selected haplos to fasta file
#################################################
filter.haplos <- function(in.flnm,out.flnm,pcnt=0.5,mnrd=NULL)
{ mnr <- 1
  if(!is.null(mnrd)) mnr <- mnrd
  lst <- read.ampl.seqs(in.flnm,mnr)
  tnr <- sum( lst$ID$nseqs)
  fl <- rep(TRUE,length(lst$seqs))
  if(!is.null(pcnt))
  {  mnrd <- round(lst$tnr*pcnt/100)
     fl <- lst$IDs$nseqs >= mnrd
  }
  if(sum(fl)==0) return()
  apc <- round(lst$ID$nseqs[fl]/sum(lst$ID$nseqs[fl])*100,2)
  nms <- paste(lst$IDs$ID[fl],lst$ID$nseqs[fl],apc,sep="|")
  seqs <- lst$seqs[fl];  
  names(seqs) <- nms
  writeXStringSet(DNAStringSet(seqs),out.flnm)
  list(df=data.frame(ID=lst$IDs$ID[fl],reads=lst$ID$nseqs[fl],Pctg=apc),
       seqs=seqs,tnr=tnr)
}


###  Lectura de la estructura de las descripciones de las muestras
###################################################################
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                    stringsAsFactors=F)
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                    stringsAsFactors=F)
RefSeqs <- as.character(readDNAStringSet(file.path(dataDir,"AmpliconRefSeqs.fna")))
idx <- sapply(samples$RefSeq.ID,function(pr) which(names(RefSeqs)==pr)[1])
samples$Primer.idx <- idx

##  Leer table de ficheros a tratar
####################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
FlTbl <- FlTbl[ FlTbl$Str=="fw", ]
nms.fw <- FlTbl$File.Name

###  Data
in.files <- FlTbl$File.Name
in.files <- sub(".PrFW.",".",nms.fw)
in.files <- file.path(joinDir,in.files)
out.files <- sub(".PrFW.",".CH01.",nms.fw)
out.files <- file.path(nt01Dir,out.files)

sink(file=file.path(repDir,"FltrConsHaplos.10-rprt.txt"))

###  Filter and save
cat("\n   Filtering consensus haplotypes by abundance")
cat("\n=================================================\n")
cat("\n  Min percentage:",ab.thr,"   min reads:",min.reads)
cat("\n\nContents by output file:\n\n")
cnr <- integer(length(in.files))
fnr <- integer(length(in.files))
nh <- integer(length(in.files))
for(i in 1:length(in.files))
{ if( !file.exists(in.files[i]) ) next
  print(out.files[i])
  lr <- filter.haplos(in.files[i],out.files[i],ab.thr,min.reads)
  df <- lr$df
  cat("\nHaplotypes:",nrow(df),"   Reads:",sum(df$reads),"\n\n")
  print(df)
  nh[i] <- nrow(df)
  cnr[i] <- lr$tnr
  fnr[i] <- sum(df$reads)

  ##  Report point mutations
  if(length(lr$seqs)>1)
  { pr.idx <- FlTbl$Pr.ID[i]
    off <- primers$FW.tpos[pr.idx]-1
    muts <- SummaryMuts.w(lr$seqs,lr$df$reads,off,RefSeqs[pr.idx][[1]])
    cat("\nObserved point mutations:\n")
    print(muts)
  } else if(length(lr$seqs) == 1)
  {
    pr.idx <- FlTbl$Pr.ID[i]
    off <- primers$FW.tpos[pr.idx]-1
    muts <- mutations(RefSeqs[samples$Primer.idx][i][[1]], lr$seqs, off, lr$df$reads)
    if (nrow(muts) != 0) {
      cat("\nObserved point mutations:\n")
      print(muts)
    }
  }
  cat("\n--------------------------------------------------\n\n")
}
fdf <- data.frame(FlTbl[,c(2,3,5)],c.reads=cnr,f.reads=fnr,
                  f.haplo=nh,pct=round(fnr/cnr*100,2))
fdf <- fdf[cnr>0,]
cat("Cutting at ",ab.thr,"%\n\n",sep="")
cat("Final yield by sample:\n\n")
print(fdf)
cat("\n\nGlobally:\n\n")
tfdf <- as.vector(apply(fdf[,4:5],2,sum))
gbl.res <- data.frame(iReads=tfdf[1],fReads=tfdf[2],
             pct2=round(tfdf[2]/tfdf[1]*100,2))
print(gbl.res)
cat("\n==================================================\n")
sink()


txt.flnm <- "FltrConsHaplos.10-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm))

###  Filter and save
cat("\n   Filtering consensus haplotypes by abundance")
cat("\n==================================================\n")
cat("\n  Min percentage:",ab.thr,"   min reads:",min.reads,"\n\n")
cat("Final yield by sample:\n\n")
print(fdf)
cat("\n\nGlobally:\n\n")
tfdf <- as.vector(apply(fdf[,4:5],2,sum))
gbl.res <- data.frame(iReads=tfdf[1],fReads=tfdf[2],
             pct2=round(tfdf[2]/tfdf[1]*100,2))
print(gbl.res)
cat("\n==================================================\n")
sink()

file.copy( file.path(repDir,txt.flnm), file.path(exportDir,txt.flnm),
           overwrite=TRUE )

