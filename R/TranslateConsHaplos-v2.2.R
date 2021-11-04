#######################################################
###     TRADUCCIO A AMINO ACIDS DELS HAPLOTIPS      ###
#######################################################

source(file.path(codeDir,"seqanalfns.aa.v4.3.R"))
source(file.path(codeDir,"seqanalfns.v4.5.R"))

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
###   la màster i de la seva frequencia poblacional, i salva les 
###   sequencies en un fitxer fasta.
########################################################################
SaveAaHaplotypes <- function(flnm,bseqs,nr)
{
  ##  Determinar dominant
  imstr <- which.max(nr)
  ##  Determinar diferències respecte la dominant
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[imstr])
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
  nms <- paste("Hpl",nm,zeroFillInt2Char(isq,4),sep=".")
  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")
  ##  Salvar a fasta
  write.fasta(AAStringSet(bseqs),flnm)
  list(bseqs=bseqs,nr=nr,nm=nm)
}


###  Recolapsar reads
Recollapse <- function(bseqs,nr)
{
  for(i in 1:(length(bseqs)-1))
  { if(nr[i]==0) next
    for(j in (i+1):length(bseqs))
    { if(nr[i] == 0 | nr[j] == 0) next
      if(bseqs[i]==bseqs[j]) 
      { nr[i] <- nr[i]+nr[j]
        nr[j] <- 0 
      } 
    }
  }
  flags <- nr > 0
  list(bseqs=bseqs[flags],nr=nr[flags])
}


###  Variabilitat en els haplotips
aa.AmplStats <- function(seqs,nr)
{
  smres <- data.frame( Hapl=integer(1),Eta=integer(1),
                       S=integer(1),Mf=numeric(1),       #Sn=numeric(1),
                       Pi=numeric(1) )
  smres$Hapl <- length(seqs)

  ###  Compute the aa frequencies and the mutations table
  pos.tbl <- aa.SeqsTbl.w(seqs,nr)
  mut.tbl <- aa.MutsTbl(pos.tbl)

  ### Total number of mutations (Eta)
  smres$Eta <- aa.TotalMutations(mut.tbl)
  cat("Total number of mutations (Eta): ",smres$Eta,"\n")

  ### Number of polymorphic sites (segregating sites, S)
  smres$S <- aa.PolymorphicSites(mut.tbl)
  cat("Number of polymorphic sites (S) : ",smres$S,"\n")

  ### Mutations frequency
  smres$Mf <- aa.TotalMutations.w(mut.tbl) / 
                (sum(nr)*nchar(seqs[1]))
  cat("Mutation frequency (Mf) : ",smres$Mf,"\n")

  ### Normalized quasispecies Shannon entropy
  # smres$Sn <- aa.QSEntropy(nr)/log(smres$Hapl)
  # cat("Normalized quasispecies Shannon entropy (Sn) : ",
       # smres$Sn,"\n")

  ###  Pairwise differences
  d.ij <- aa.PairwiseDiffs.w(seqs,nr)

  ###  Mean number of amino acid differences
  cat("Mean number of pairwise amino acid differences (k) : ",
       MeanAaDiffs.w(d.ij,nr),"\n")

  ###  Amino acid diversity (Pi)
  smres$Pi <- MeanAaDiffs.w(d.ij,nr)/nchar(seqs[1])
  cat("Amino acid diversity (Pi) : ",smres$Pi,"\n")

  smres
}


###   MAIN LINE
#######################################################################

library(Biostrings)
library(RColorBrewer)
cls <- brewer.pal(8,"Dark2")

###  Llegim l'estructura de descripció de mostres
##################################################
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                    stringsAsFactors=F)

###  Llegim fitxers globals de HCV
##################################################
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                    stringsAsFactors=F)
###  Taula d'ORFs
orfs <- data.frame(nt.start=primers$FW.tpos,nt.end=primers$RV.tpos,
                   aa.ipos=primers$Aa.ipos,aa.lpos=primers$Aa.lpos)

###  Localitzar amplicons de l'experiment
##################################################
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
if( any(is.na(idx)) )
  stop("No primers descriptor found for samples: ",
        paste(which(is.na(idx)),collapse=", "))
samples$Primer.idx <- idx

##  Llegim la taula de fitxers a tractar
##########################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))


###  Generar noms de fasta, d'entrada i de sortida
FlTbl <- FlTbl[FlTbl$Str=="fw",]
oflnms <- FlTbl$File.Name

flnms <- sub(".PrFW.",".CH05.",oflnms)
flnms <- file.path(ntDir,flnms)

out.flnms <- sub(".PrFW.",".AaCH05.",oflnms)
out.flnms <- file.path(aaDir,out.flnms)

###  Inici d'informe
######################
sink(file=file.path(repDir,"HaplosTranslation-Rprt.txt"))

cat("\n   HAPLOTYPES TRANSLATION")
cat("\n============================\n")

n <- length(flnms)
rlst <- vector("list",n)
rdf <- data.frame(reads=integer(n),nt.haplos=integer(n),aa.haplos=integer(n))

res <- data.frame( Hapl=integer(n),Eta=integer(n),S=integer(n),Mf=numeric(n),
                   Pi=numeric(n) )
fmast <- numeric(n)
for(i in 1:n)
{  pr.idx <- FlTbl$Pr.ID[i]

   if( ! file.exists(flnms[i]) ) next

   cat("\nProcessing file:")
   cat("\n   ",flnms[i])
   cat("\n    nt positions ",orfs$nt.start[pr.idx]," to ",
       orfs$nt.end[pr.idx],sep="")
   cat("\n    aa positions ",orfs$aa.ipos[pr.idx]," to ",
       orfs$aa.lpos[pr.idx],"\n",sep="")
   off <- orfs$aa.ipos[pr.idx]-1
   res$Hapl[i] <- 1

   ###  Read nucleotide haplotypes
   lst <- read.ampl.seqs(flnms[i],mnr=1)
   nms <- names(lst$seqs)
   rdf[i,1] <- sum(lst$IDs$nseqs)
   rdf[i,2] <- length(lst$seqs)
   cat("\n    Nucleotide haplotypes:",rdf[i,2])
   cat("\n    Nucleotide reads:",rdf[i,1],"\n")

   ###  Trim to orf
   i0 <- orfs$nt.start[pr.idx]-primers$FW.tpos[pr.idx]+1
   i1 <- i0+(orfs$nt.end[pr.idx]-orfs$nt.start[pr.idx])
   if(i0<1 | i1>nchar(lst$seqs[1]))
      stop("Wrong orf coordinates.\n")
   if( i0 > 1 || i1 < nchar(lst$seqs[1]))
      lst$seqs <- substring(lst$seqs,i0,i1) 

   ###  Translate to aa
   aa.seq <- translate(DNAStringSet(lst$seqs))
   aa.seq <- as.character(aa.seq)
   names(aa.seq) <- nms

   ###  Recollapse
   lst <- list(bseqs=aa.seq,nr=lst$IDs$nseqs)
   if(length(aa.seq)>1)
     lst <- Recollapse(aa.seq,lst$nr)

   ###  Salvar fasta d'amino acids
   lst <- SaveAaHaplotypes(out.flnms[i],lst$bseqs,lst$nr)

   ###  Informar del resultat de la traducció i posterior recollapsat
   nr <- lst$nr
   cat("\nTable of aa haplotypes differences:")
   print(table(lst$nm))
   cat("\nTable of reads differences:\n")
   print(tapply(nr,as.factor(lst$nm),sum))
   rdf[i,3] <- length(nr)
   cat("\nFinal recolapsed aa haplotypes:",length(nr))
   cat("\nTotal aa reads:",sum(nr),"\n")

   ###  Taula de llocs polimòrfics
   if(length(nr)>1)
   { cat("\nTable of aa polymorphic sites:\n")
     sm <- aa.SummaryMuts.w(lst$bseqs,lst$nr,off)  
     rlst[[i]] <- data.frame(sm)
     cat(" pos: list(aa counts,...)\n")
     print.aa.polym(sm[,(1:(naas+1))])
     cat("\n pos: list(aa frequency%,...)\n")
     print.aa.polym(sm[,c(1,(1:naas)+naas+1)])
     cat("\n")
     ###  Variabilitat en els haplotips
     res[i,] <- aa.AmplStats(lst$bseqs,lst$nr)
     fmast[i] <- round(lst$nr[1]/sum(lst$nr)*100,2)
   } else 
       cat("No aa mutations observed\n")

   cat("\n-------------------------------------------------------\n")
}

sink()


txt.flnm <- "HaplosTranslation-SummaryRprt.txt"
sink(file=file.path(repDir,txt.flnm),split=TRUE)
cat("\n   HAPLOTYPES TRANSLATION")
cat("\n============================\n")
cat("\nGLOBAL SUMMARY\n\n")
print(rdf)
cat("\n")
df <- data.frame(Pat.ID=FlTbl$Pat.ID,Ampl.Nm=FlTbl$Ampl.Nm,
                 signif(res,4),Mstr=fmast)
print(df)
sink()

file.copy( file.path(repDir,txt.flnm), file.path(exportDir,txt.flnm),
           overwrite=TRUE )
