#############################################################
###    NUCLEOTIDE SEQUENCE ANALYSIS, UTILITY FUNCTIONS    ###
###    - Pairwise differences excluding gaps              ###
###                                                       ###
###             source("seqanalfns.v4.0.txt")             ###
#############################################################
library(Biostrings)
library(RColorBrewer)
cls <- brewer.pal(8,"Dark2")
nuc.nms <- c("A","C","G","T","Gap")


###  Table of nucleotide frequencies at position i
PosTbl <- function(i,seqs,w)
{ res <- c(A=0,C=0,G=0,T=0,"-"=0)
  tbl <- table(substr(seqs,i,i))
  res[names(tbl)] <- tbl
  res
}


###  Table of population nucleotide frequencies at position i
###    given sequence frequencies
PosTbl.w <- function(i,seqs,w)
{ alph <- c("A","C","G","T","-")
  res <- c(A=0,C=0,G=0,T=0,Gap=0)
  lets <- substr(seqs,i,i)
  for(i in 1:5)
    res[i] <- sum(ifelse(lets==alph[i],1,0)*w)
  res
}


###  Matrix of nucleotide frequencies per position
SeqsTbl <- function(seqs)
{ tbl <- t(sapply(1:nchar(seqs[1]),function(i) PosTbl(i,seqs)))
  rownames(tbl) <- 1:nchar(seqs[1])
  tbl
}

###  Matrix of population nucleotide frequencies per position
###    given sequence frequencies
SeqsTbl.w <- function(seqs,w)
{ tbl <- t(sapply(1:nchar(seqs[1]),function(i) PosTbl.w(i,seqs,w)))
  rownames(tbl) <- 1:nchar(seqs[1])
  tbl
}


###  Matrix of mutations frequencies per position
###  (the most frequent is set to 0 leaving only the mutations)
MutsTbl <- function(seq.tbl)
{ j <- apply(seq.tbl,1,which.max)
  seq.tbl[cbind(1:nrow(seq.tbl),j)] <- 0
  seq.tbl
}


###  Consensus sequence build from the matrix of nucleotide 
###    frequencies per position. The most frequent nucleotide 
###    at each position is taken.
ConsSeq <- function(seq.tbl)
{ j <- apply(seq.tbl,1,which.max)
  paste(colnames(seq.tbl)[j],collapse="")
}


###  Matrix of summary mutations
SummaryMuts <- function(seqs,off=0)
{ pos.tbl <- SeqsTbl(seqs)
  mut.tbl <- MutsTbl(pos.tbl)
  flags <- apply(mut.tbl,1,sum)>0
  pos <- which(flags)
  res <- cbind(pos=pos+off,pos.tbl[flags,],
               round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
  colnames(res) <- c("pos",paste("f",nuc.nms,sep= ""),
                     paste("p",nuc.nms,sep= ""))
  rownames(res) <- 1:nrow(res)
  res
}


###  Matrix of summary population mutations
###    given sequence frequencies
SummaryMuts.w <- function(seqs,w,off=0)
{ if(length(seqs)<2)
    return(NULL)
  pos.tbl <- SeqsTbl.w(seqs,w)
  mut.tbl <- MutsTbl(pos.tbl)
  flags <- apply(mut.tbl,1,sum)>0
  pos <- which(flags)
  if(length(pos)==1)
  { res <- c(pos=pos+off,pos.tbl[flags,],
               round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
    res <- matrix(res,nrow=1)
  } else {	  
    res <- cbind(pos=pos+off,pos.tbl[flags,],
                 round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
  }
  rownames(res) <- 1:nrow(res)
  colnames(res) <- c("pos",paste("f",nuc.nms,sep= ""),
                     paste("p",nuc.nms,sep= ""))
  data.frame(res)
}


###  Segregating sites: Number of sites with mutations
###    (independent of sequence weights)
SegSites <- function(seqs)
{ nrow(SummaryMuts(seqs)) }


###  Segregating sites: Number of sites with mutations
###    given sequence frequencies
SegSites.w <- function(seqs,w)
{ nrow(SummaryMuts.w(seqs)) }


### Total number of mutations (Eta)
###  (independent of sequence weights)
TotalMutations <- function(mut.tbl)
  sum(apply(mut.tbl,1,function(x) sum(x>0)))


### Total number of population mutations
TotalMutations.w <- function(mut.tbl,w)
  sum( apply(mut.tbl,1,function(x) sum(x)) )


### Number of polymorphic sites (segregating sites) (S)
###  (independent of sequence weights)
PolymorphicSites <- function(mut.tbl)
  sum( apply(mut.tbl,1,sum) > 0 )


###  Nucleotide differences between two sequences (excluding gaps) 
PairDiffs <- function(seq1,seq2)
{ s1 <- strsplit(seq1,split="")[[1]]
  s2 <- strsplit(seq2,split="")[[1]]
  sum( s1 != s2 ) - 
     sum( ((s1=="-")&(s2!="-")) | ((s1!="-")&(s2=="-")) )
}


###  Matrix of pairwise differences
PairwiseDiffs <- function(seqs)
{ n <- length(seqs)
  d <- matrix(0,n,n)
  if(length(seqs)<2) return(d)
  for(i in 1:(n-1))
    for(j in (i+1):n)
       d[i,j] <- PairDiffs(seqs[i],seqs[j])
  d+t(d)
}


###  Matrix of population pairwise differences
###    given sequence frequencies
PairwiseDiffs.w <- function(seqs,w)
{ n <- length(seqs)
  d <- matrix(0,n,n)
  if(length(seqs)<2) return(d)
  for(i in 1:(n-1))
    for(j in (i+1):n)
       d[i,j] <- PairDiffs(seqs[i],seqs[j])*w[i]*w[j]
  d+t(d)
}


###  Mean number of nucleotide differences
MeanNuclDiffs <- function(d)
{ n <- nrow(d)
  if(n<2) return(0)
  k <- 0
  for(i in 1:(n-1))
    for(j in (i+1):n)
       k <- k + d[i,j]
  k*2/(n*(n-1))
}
 
 
###  Mean number of population nucleotide differences
###    given sequence frequencies
MeanNuclDiffs.w <- function(d,w)
{ n <- nrow(d)
  if(n<2) return(0)
  nseqs <- sum(w)
  k <- 0
  for(i in 1:(n-1))
    for(j in (i+1):n)
       k <- k + d[i,j]
  k*2/(nseqs*(nseqs-1))
}


###  Information content at a site
###    given the nucleotides frequency
InfContent <- function(v)
{ v <- v/sum(v)
  lgv <- ifelse(v==0,0,log2(v))
  2+sum(v*lgv)
}


###  Information content at a site given the absolute mutations
###    frequency, and the total nunber of population sequences
InfContent.mut <- function(v,tseq)
{ v <- c(v/tseq,(tseq-sum(v))/tseq)
  lgv <- ifelse(v==0,0,log2(v))
  2+sum(v*lgv)
}


###  Information content at a site
###    given the mutations percentages
InfContent.pmut <- function(v)
{ v <- v/100
  v <- c(v,1-sum(v))
  lgv <- ifelse(v==0,0,log2(v))
  2+sum(v*lgv)
}


###  Quasispecies Shannon entropy
QSEntropy <- function(w)
{ w <- w/sum(w)
  lgw <- ifelse(w==0,0,log(w))
  -sum(w*lgw)
}


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


### After right and/or left trimming recollapse identical sequences
##################################################################### 
CollapseSeqs <- function(seqs,IDs)
{ if(length(seqs)<2)
    return( list(seqs=seqs,IDs=IDs) )
  for(i in 1:(length(seqs)-1))
  { if(IDs$nseqs[i]==0) next
    for(j in (i+1):length(seqs))
    { if(IDs$nseqs[i] == 0 | IDs$nseqs[j] == 0) next
      if(seqs[i]==seqs[j]) 
      { IDs$nseqs[i] <- IDs$nseqs[i]+IDs$nseqs[j]
        IDs$nseqs[j] <- 0 
      } 
    }
  }
  flags <- IDs$nseqs > 0
  list(seqs=seqs[flags],IDs=IDs[flags,])
}

             ###   END OF seqanalfns.v4.5.R   ###
 