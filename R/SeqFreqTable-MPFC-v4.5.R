##############################################################################
###                 MODUL: SeqFreqTable-FullCorrect          
###
###  Parteix dels fasta de reads obtesos del MiSeq un cop demultiplexats 
###    i retallats.
###
###  Es corregeixen gaps i Ns per la RefSeq (o la màster)
###  Els haplotipus amb més de max.gaps o max.Ns s’eliminen
###  Si max.gaps es 0 s’eliminen tots els reads amb gaps
###  Si max.Ns es 0 s’eliminen tots els reads amb Ns
###
###  Es genera un informe sobre els resultats obtesos.
##############################################################################

source("./R/seqanalfns.v4.5.R")
library(data.table)

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


###  Funcions de diferencies entre seqüències aliniades
########################################################

##  Diferències posició a posició
n.difs <- function(seq1,seq2)
{ sum( strsplit(seq1,split="")[[1]] != strsplit(seq2,split="")[[1]] ) }

##  Gaps a un i altre costat
n.gaps <- function(seq1,seq2)
{ seq1 <- strsplit(seq1,split="")[[1]]
seq2 <- strsplit(seq2,split="")[[1]]
sum( (seq1=="-" & seq2!="-") | (seq1=="-" & seq2!="-") )
}

##  Diferències menys gaps
n.edit <- function(seq1,seq2)
{ seq1 <- strsplit(seq1,split="")[[1]]
seq2 <- strsplit(seq2,split="")[[1]]
sum(seq1!=seq2)-sum((seq1=="-" & seq2!="-")|(seq1=="-" & seq2!="-"))
}

##  Identitat % entre dues seqüències
pid2 <- function(seq1,seq2)
  100*(1-n.edit(seq1,seq2)/nchar(seq1))



###  Compta les Ns de cada sequencia
######################################
has.Ns <- function(seq)
  as.vector(apply(as.matrix(seq),1,function(ch) sum(ch=="N")))


###  Determina coincidències d'un segment dins d'un conjunt de
###    fragmnents de seqüències.
#################################################################
procBySubstr <- function(seqs,patt,mm.mx,ambig)
{
  ##  Buscar coincidències. No indels. Ambigüetats depenent de 'ambig'
  matches <- vmatchPattern(patt,subject=seqs,fixed=!ambig,max.mismatch=mm.mx,
                           with.indels=TRUE) ## CAMBIO. with.indels = TRUE
  ##  elementLengths() és de dimensió completa amb les posicions de
  ##    no coincidència a 0. Evitar coincidències múltiples.
  flags <- elementLengths(matches)==1
  ##  Aquest pas és delicat
  pos <- unlist(startIndex(matches)[flags])
  end <- unlist(endIndex(matches)[flags])
  list(flags=flags,pos=pos,end=end)
}


###  Table of population nucleotide frequencies at position i
###    given sequence frequencies
###############################################################
PosTbl.w <- function(i,seqs,w)
{ alph <- c("A","C","G","T")
res <- c(A=0,C=0,G=0,T=0)
lets <- substr(seqs,i,i)
for(i in 1:4)
  res[i] <- sum(ifelse(lets==alph[i],1,0)*w)
res
}


###  Matrix of population nucleotide frequencies per position
###    given sequence frequencies
##################################################################
SeqsTbl.w <- function(seqs,w)
{ tbl <- t(sapply(1:nchar(seqs[1]),function(i) PosTbl.w(i,seqs,w)))
rownames(tbl) <- 1:nchar(seqs[1])
tbl
}


###  Consensus sequence build from the matrix of nucleotide
###    frequencies per position. The most frequent nucleotide
###    at each position is taken.
##################################################################
ConsSeq <- function(seq.tbl)
{ j <- apply(seq.tbl,1,which.max)
paste(colnames(seq.tbl)[j],collapse="")
}


###  Recolapsar reads
######################
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


###  Corregir les Ns de sequencies amb menys de max.Ns
###   Si hi ha més de max.Ns en una sequencia la deixa inalterada.
####################################################################
CorrNs <- function(cseqs, RefSeq)
{
  seqlist <- list()
  nslist <- list()
  for (i in 1:length(as.matrix(cseqs))) {
    ## Sequencies a corregir (fins a un màxim de max.Ns)
    ns <- has.Ns(DNAStringSet(cseqs[i]))
    flc <- ns>0 & ns<=max.Ns
    if (sum(flc)==0) {
      seqlist[[i]] <- as.character(cseqs[i])
      nslist[[i]] <- as.integer(ns)
    } else {
      ## Posició a corregir
      sqs <- DNAStringSet(cseqs[i][flc])
      cm <- consensusMatrix(sqs)
      pos <- which(cm["N",]>0)
      
      ## Base en aquesta posició per la RefSeq
      let <- strsplit(RefSeq[i], split="")[[1]][pos]
      for (j in 1:length(pos)) {
        ## vector de lògics indicant les seqs que tenen una N en aquesta pos
        fl <- letter(as.character(sqs), pos[j]) == "N"
        ## matriu de lògics indicant les posicions on hi ha Ns a corregir
        at <- matrix(FALSE, nrow=sum(fl), ncol=nchar(sqs[1]))
        at[cbind(1:nrow(at), pos[j])] <- TRUE
        ## correcio
        sqs[fl] <- replaceLetterAt(DNAStringSet(sqs[fl]), at=at, letter=rep(let[j], sum(fl)))
      }
      seqlist[[i]] <- as.character(sqs)
      nslist[[i]] <- as.integer(ns)
    }
  }
  list(cseqs=as.character(seqlist),ns=as.integer(nslist))
}

###  Corregeix tots els gaps pel contingut de la sequencia de referencia
###   Si hi ha més de max.gaps en una sequencia la deixa inalterada.
###########################################################################
CorrectGaps <- function(bseqs,ref.seq)
{ rf.lt <- strsplit(ref.seq,split="")[[1]]
ngap <- alphabetFrequency(DNAStringSet(bseqs))[,"-"]
for(i in 1:length(bseqs))
{ if(ngap[i]==0 | ngap[i]>max.gaps) next
  lt <- strsplit(as.character(bseqs[i]),split="")[[1]]
  idx <- which(lt=="-")
  lt[idx] <- rf.lt[idx]
  bseqs[i] <- paste(lt,collapse="")
}
return(list(bseqs=bseqs,ngap=ngap))
}


###  Dóna nom als haplotipus en funció del nombre de mutacions respecte
###   la referència (que pot ser la màster) i de la seva frequencia
###   poblacional, i salva les sequencies en un fitxer fasta.
########################################################################
SaveHaplotypes <- function(ifl,bseqs,nr)
{
  ##  Nom genèric haplotipus
  code <- paste(FlTbl$Pat.ID[ifl],FlTbl$Ampl.Nm[ifl],FlTbl$Str[ifl],sep=".")
  
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
  flnm <- sub("TCA.454Reads","nt.FCHaplos",FlTbl$File.Name[ifl])
  write.fasta(DNAStringSet(bseqs),file.path(filtDir,flnm))
  
  list(bseqs=bseqs,nr=nr,nm=nm)
}


###  Determinar màster i eventualment corregir-la per la RefSeq.
####################################################################
CorrectMaster <- function(bseqs,nr,RefSeq=NULL)
{
  ##  Si no hi ha reads
  if(!length(bseqs))
  { cat("\nNo reprentative sequences found.\n")
    cat("\n-------------------------------------------------\n\n")
    return(NULL)  ## Conjunt buit.
  }
  
  ##  Si no hi ha sequencia de referencia tornar
  if(is.null(RefSeq)) 
    return( list(bseqs=bseqs, nr=nr, mseq=bseqs[1]) )
  
  ##  Alinear màster respecte la seq de referència
  psa <- pairwiseAlignment(pattern=bseqs[1],subject=RefSeq)
  mseq <- as.character(aligned(psa))[1]
  
  ###  Corregir Ns i gaps en la màster per la RefSeq
  lt1 <- strsplit(mseq,split="")[[1]]
  idx <- c(which(lt1=="-"),which(lt1=="N"))
  if(length(idx))
  { cat("\nPositions corrected on the master:\n")
    print(idx)
    lt2 <- strsplit(RefSeq,split="")[[1]]
    lt1[idx] <- lt2[idx]
    bseqs[1] <- paste(lt1,collapse="")
  }
  list(bseqs=bseqs, nr=nr, mseq=bseqs[1])
}


###  Alinea haplos respecte la màster (corregida) i els filtra i
###   corregeix segons paràmetres del programa
####################################################################
FilterCorrectHaplos <- function(ifl,bseqs,nr,mnrd=2,RefSeq=NULL)
{
  low.freq <- 0
  res <- data.frame(tmdifs=integer(1),tmNs=integer(1),tmdels=integer(1),
                    corr=integer(1),f.reads=integer(1),f.haplo=integer(1))
  # res <- data.frame(Ns=integer(1),tmdifs=integer(1),tmdels=integer(1),
  # dels=integer(1),ins=integer(1),corr=integer(1),
  # f.reads=integer(1),f.haplo=integer(1))
  cons.seq <- ""
  names(cons.seq) <- paste(FlTbl$Pat.ID[ifl],FlTbl$Ampl.Nm[ifl],FlTbl$Str[ifl],
                           sep=".")
  
  ##  Alinear respecte la RefSeq o la master (prèviament corregida)
  mseq <- RefSeq
  if(is.null(mseq)) mseq <- bseqs[1]  ##  usar la master
  psa <- pairwiseAlignment(pattern=bseqs,subject=mseq)
  
  
  ##  Invocar aligned() és fonamental: torna les sequencies amb gaps on
  ##   hi ha insercions respecte la ref, i amb les posicions que donen
  ##   gaps en la ref eliminades, el que comportarà gaps a final de cadena.
  bseqs <- as.character(pattern(psa))
  refs_aln <- as.character(subject(psa))
  
  ##  Informar de sequencies amb massa differències
  ndiffs <- nmismatch(psa)
  dif.fl <- ndiffs <= max.diffs         
  htmd <- sum(!dif.fl)
  res$tmdifs <- rtmd <- sum(nr[!dif.fl])    ###
  if(htmd)
  { cat("\nExcluded haplotypes with too many differences:",sum(!dif.fl))
    cat("  (",rtmd," reads)",sep="")
  }
  
  ##  Informar de delecions (gaps en els reads)
  dels <- deletion(psa)
  del.fl <- sapply(dels,function(l) length(l)>0)
  if(sum(del.fl))
  { cat("\nHaplotypes with deletions:",sum(del.fl) )
    cat("  (",sum(nr[del.fl])," reads)",sep="")
  }
  
  ##  Informar d'insercions (gaps en la referencia, pos eliminades en reads)
  ins <- insertion(psa)
  ins.fl <- sapply(ins,function(l) length(l)>0)
  if(sum(ins.fl))
  { cat("\nHaplotypes with insertions:",sum(ins.fl) )
    cat("  (",sum(nr[ins.fl])," reads)",sep="")
  }
  
  ##  Corregir Ns
  cNs <- 0
  lst <- CorrNs(bseqs,refs_aln)
  bseqs <- lst$cseqs
  fl.Ns.corr <- lst$ns>0 & lst$ns<=max.Ns
  fl.tmn <- lst$ns>max.Ns
  htmn <- sum(fl.tmn)
  res$tmNs <- rtmn <- sum(nr[fl.tmn])       ###
  if(length(bseqs) & sum(lst$ns))
  { cNs <- sum(fl.Ns.corr) # nombre d’haplotipus corregit
  res$corr <- sum(nr[fl.Ns.corr])
  cat("\n\nTable of Ns occurrences by haplotype:")
  print(table(lst$ns))
  cat("Haplotypes with Ns corrected:",cNs)
  cat("  (",sum(nr[fl.Ns.corr])," reads)",sep="")
  cat("\nTotal Ns in haplotypes corrected:",sum(lst$ns[fl.Ns.corr]))
  cat("\nTotal Ns in reads corrected:",
      sum(lst$ns[fl.Ns.corr]*nr[fl.Ns.corr]))
  if(htmn)
    cat("\nExcluded ",htmn," haplotypes (",rtmn," reads) with too many Ns.",
        sep="")
  }
  
  ##  Corregir gaps (provinents de deleccions)
  fl.tmg <- rep(FALSE,length(bseqs))
  if(sum(del.fl)+sum(ins.fl))
  { 
    ngap <- alphabetFrequency(DNAStringSet(bseqs))[,"-"]
  #bseqs <- lst$bseqs
  #ngap <- lst$ngap
    fl.ng <- ngap > 0 & ngap <= max.gaps
    fl.tmg <- ngap > max.gaps
  #htmg <- sum(fl.tmg)                 
  #res$tmdels <- rtmg <- sum(nr[fl.tmg])   ###
  #cat("\n\nTable of gap occurrences by haplotype:")
  #print(table(lst$ngap))
  #cat("\nHaplotypes with gaps corrected:",sum(fl.ng))
  #cat("  (",sum(nr[fl.ng])," reads)",sep="")
  #cat("\nTotal gaps in haplotypes corrected:",sum(ngap[fl.ng]))
  #cat("\nTotal gaps in reads corrected:",sum(ngap[fl.ng]*nr[fl.ng]))
  #if(htmg)
  #cat("\nExcluded ",htmg," haplotypes (",rtmg,
  #      " reads) with too many gaps.",sep="")
  #res$corr <- sum(nr[fl.Ns.corr | fl.ng])
  }
  
  cat("\n\nTotal haplotypes excluded by Ns, gaps or differences:",
      sum( fl.tmn | fl.tmg | !dif.fl ) )
  cat("  (",sum(nr[fl.tmn | fl.tmg | !dif.fl])," reads)",sep="")
  
  ##  Eliminar els haplotipus que no passen el filtre
  bseqs <- bseqs[!fl.tmn & !fl.tmg & dif.fl]
  nr <- nr[!fl.tmn & !fl.tmg & dif.fl]
  
  ##  Si no hi ha sequencies vàlides
  if(length(bseqs)==0)
  { cat("\nNo sequences passed the filter of Ns + indels + max diffs.\n")
    cat("\n----------------------------------------------------------\n\n")
    return(list(res=res,mseq=NULL,nrmseq=0,low.freq=low.freq))  ## Buit.
  }
  
  if(length(nr)>1)
  {
    ##  Recolapsar
    lst <- Recollapse(bseqs,nr)
    bseqs <- lst$bseqs
    nr <- lst$nr
    ##  AQUÍ LA MASTER POT SER UNA ALTRA DESPRÉS DE LES REPARACIONS
    o <- order(nr,decreasing=TRUE)
    bseqs <- bseqs[o]
    nr <- nr[o]
    if( mseq != bseqs[1] )
    { cat("\n  Warning: a new master emerged after the repairs")
      cat("\n           with frequency",nr[1],"reads")
      mseq <- bseqs[1]
    }
    
    ##  Post-filtre per nombre de reads
    fl <-  nr>=mnrd
    if(sum(!fl))
    { cat("\nHaplotypes with less than",mnrd,"reads, excluded:",sum(!fl))
      cat("\n representing",sum(nr[!fl]),"reads.")
      low.freq <- sum(!fl)                  ###
    }
    nr <- nr[fl]
    bseqs <- bseqs[fl]
    if(length(bseqs)==0)
    { cat("\nNo haplotypes found with at least",min.reads,"reads.\n")
      cat("\n-------------------------------------------------\n\n")
      return(list(res=res,mseq=NULL,nrmseq=0,low.freq=low.freq))  ## Buit.
    }
    
    ##  Generar la consens poblacional.
    consMat <- SeqsTbl.w(bseqs,nr)
    cons.seq <- ConsSeq(consMat)
    names(cons.seq) <- "000|0"
    
    ##  Diferències entre la master i la consens  
    cmdifs <- n.edit(bseqs[1],cons.seq)
    cat("\n\nConsensus – master differences:",cmdifs)
  }
  cat("\nFinal reads for the master sequence:",nr[1])
  
  ## Ordenar per mutacions i freqüències, anomenar, i salvar haplotips
  lst <- SaveHaplotypes(ifl,bseqs,nr)
  nr <- lst$nr
  cat("\n\nTable of haplotypes differences:")
  print(table(lst$nm))
  cat("\nTable of reads differences:\n")
  print(tapply(nr,as.factor(lst$nm),sum))
  
  cat("\nFinal recolapsed haplotypes passing the filters:",length(bseqs))
  cat("\nFinal reads for the master sequence:",nr[1])
  cat("\nFinal reads passing the filters:",sum(nr))
  cat("\n-------------------------------------------------\n\n")
  
  res$f.reads <- sum(nr)
  res$f.haplo <- length(bseqs)
  list(res=res,mseq=bseqs[1],nrmseq=nr[1],low.freq=low.freq)
}


############################################################################

###  MAIN LINE  ###

cls <- brewer.pal(8,"Dark2")

###  Llegim l'estructura de descripció de mostres
##################################################
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)

###  Localitzar amplicons de l'experiment
##################################################
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
if( any(is.na(idx)) )
  stop("No primers descriptor found for samples: ",
       paste(which(is.na(idx)),collapse=", "))
samples$Primer.idx <- idx

###  Carregar sequencies de referencia generiques per cada mostra
###################################################################
RefSeqs <- as.character(
  readDNAStringSet(file.path(dataDir,"AmpliconRefSeqs.fna")))
idx <- sapply(samples$RefSeq.ID,function(pr) which(names(RefSeqs)==pr)[1])
if( any(is.na(idx)) )
  stop("No reference sequence found for samples: ",
       paste(which(is.na(idx)),collapse=", "))
###  Vincular una RefSeq a cada mostra
samples$RefSeq.idx <- idx

##  Llegim la taula de fitxers a tractar
##########################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
FlTbl$RefSeq.idx <- 0
for(i in 1:nrow(FlTbl))
{ k <- which(samples$Patient.ID == FlTbl$Pat.ID[i] &
               samples$Primer.ID == FlTbl$Ampl.Nm[i])
FlTbl$RefSeq.idx[i] <- samples$RefSeq.idx[k]
}  

##  Inici d'informe
#####################
sink(file=file.path(repDir,"FCReadsFilter-Rprt.txt"))
cat("\n   = = =   FILTERING READS FROM THE MiSeq  = = =")
cat("\n===================================================\n")
cat("\n Thresholds & conditions imposed:")
if(is.null(RefSeqs))
{ cat("\n    No reference sequence. Using master for alignment and repairs.")
} else {
  if(ref.type=="generic")
    cat("\n    Generic reference sequence for alignment and repairs.")
  if(ref.type=="consensus")
    cat("\n    Consensus reference sequence for alignment and repairs.")
}
cat("\n    Max differences tolerated:",max.diffs)
cat("\n    Ns and gaps corrected as per RefSeq contents.")
cat("\n    Max number of Ns tolerated:",max.Ns)
cat("\n    Max number of gaps tolerated:",max.gaps)
cat("\n    Post repair min reads per haplotype:",min.reads)

cat("\n\n=======================================================\n")

# res <- data.frame(i.reads=integer(),shorts=integer(),no.pr=integer(),
# low.freq=integer(),
# Ns=integer(),dels=integer(),ins=integer(),corr=integer(),
# f.reads=integer(),f.haplo=integer())

res <- data.frame(i.reads=integer(),low.freq=integer(),
                  tmdifs=integer(),tmNs=integer(),tmdels=integer(),
                  corr=integer(),f.reads=integer(),f.haplo=integer())

##  Loop sobre fitxers a analitzar
#######################################
ampl.set <- unique(sort(FlTbl$Ampl.Nm))
mseqs <- character(nrow(FlTbl))
for( i in 1:nrow(FlTbl) )
{
  ##  Sequencia de referencia pels reads nets del fitxer
  RefSeq <- NULL
  pr.id <- FlTbl$Pr.ID[i]
  if( !is.null(RefSeqs) ) 
    RefSeq <- RefSeqs[FlTbl$RefSeq.idx[i]]
  
  ##  Llegir fitxer fasta
  flnm <- FlTbl$File.Name[i]
  if( ! file.exists(file.path(trimDir,flnm)) ) next
  lst <- read.ampl.seqs(file.path(trimDir,flnm),mnr=1)
  bseqs <- lst$seqs
  nr <- lst$IDs$nseqs
  
  cat("\nFile: ",flnm)
  cat("\nPatient:",FlTbl$Pat.ID[i],"   Ampl:",FlTbl$Ampl.Nm[i],
      "   Strand:",FlTbl$Str[i])
  cat("\n\nHaplotypes ",length(nr),", reads ",sum(nr),sep="","\n")
  res[i,"i.reads"] <- sum(nr)
  
  ##  Identificar màster i corregir-la per la RefSeq si és "generic"
  rseq <- NULL
  if(!is.null(RefSeq) & ref.type=="generic") rseq <- RefSeq
  lst <- CorrectMaster(bseqs,nr,rseq)
  if(is.null(lst)) 
    next
  ##  Guardar-la
  mseqs[i] <- lst$mseq
  
  ##  Filtrar i corregir haplos per la màster si no hi ha cap RefSeq, 
  ##    per la màster ja corregida si la RefSeq és "generic", o
  ##    per la mateixa RefSeq si és "consensus"
  rseq <- NULL
  if(ref.type=="consensus") rseq <- RefSeq
  lst <- FilterCorrectHaplos(i,lst$bseqs,lst$nr,mnrd=min.reads,RefSeq=rseq)
  res[i,"low.freq"] <- lst$low.freq
  res[i,3:8] <- lst$res
  if(lst$res$f.reads > 0 & !is.null(lst$mseq))   ###  27/10/2016
    mseqs[i] <- lst$mseq
  
}      ##  Final de loop sobre amplicons

res <- res[ !is.na(res[,1]), ]

sink()

###   SUMMARY and PLOTS

txt.flnm <- "FCReadsFilter-SummaryRprt.txt"
sink(file=file.path(repDir,txt.flnm),split=TRUE)
cat("\n   = = =   FILTERING READS FROM THE 454 GS-FLX  = = =")
cat("\n========================================================\n")
cat("\n Thresholds & conditions imposed:")
if(is.null(RefSeqs))
{ cat("\n    No reference sequence. Using master for alignment and repairs.")
} else {
  if(ref.type=="generic")
    cat("\n    Generic reference sequence for alignment and repairs.")
  if(ref.type=="consensus")
    cat("\n    Consensus reference sequence for alignment and repairs.")
}
cat("\n    Max differences tolerated:",max.diffs)
cat("\n    Ns and gaps corrected as per RefSeq contents.")
cat("\n    Max number of Ns tolerated:",max.Ns)
cat("\n    Max number of gaps tolerated:",max.gaps)
cat("\n    Post repair min reads per haplotype:",min.reads)

cat("\n\n=======================================================\n")

cat("\nSummary of results:\n\n")
print(res)

t.res <- apply(res[,c(1:7)],2,sum)
cat("\n")
print(t.res)
t.res <- round(t.res/t.res[1]*100,2)
cat("\n")
print(t.res)

cat("\n\nYield in filtering reads, by file:\n\n")
FlTbl <- FlTbl[as.integer(rownames(res)),]
df <- data.frame(Pat.ID=FlTbl$Pat.ID,Ampl.Nm=FlTbl$Ampl.Nm,Str=FlTbl$Str,
                 i.Reads=res$i.reads,f.Haplos=res$f.haplo,f.Reads=res$f.reads)
print(df)
cat("\n")

df2 <- data.frame(Pat.ID=FlTbl$Pat.ID,Ampl.Nm=FlTbl$Ampl.Nm,Str=FlTbl$Str,
                  i.Reads=res$i.reads,f.Reads=res$f.reads,
                  yield=round(res$f.reads/res$i.reads*100,2),
                  pc.corr=round(res$corr/res$f.reads*100,2))
df2$pc.corr[df$yield==0] <- 0
print(df2)

##  Mostres perdudes.
nlst <- sum(df$i.Reads>0 & df$f.Reads==0)
if(nlst)
{ cat("\n\nSamples lost in the filters:  ",nlst," (",
      round(nlst/nrow(df)*100,2),"%)\n",sep="")
  fl <- df$f.Reads>0
  cat("\nGlobal table by files:\n")
  print(table(df$Ampl.Nm[fl],df$Pat.ID[fl]))
}

##  Taula de reads per amplicó/pacient (FW+RV)
tbl2 <- aggregate(df$f.Reads,list(df$Ampl.Nm,df$Pat.ID),sum)
mt <- matrix(0,nrow=length(unique(tbl2[,1])),ncol=length(unique(tbl2[,2])))
rownames(mt) <- unique(tbl2[,1])
colnames(mt) <- unique(tbl2[,2])
mt[as.matrix(tbl2[,1:2])] <- tbl2[,3]
cat("\nGlobal table by reads (FW+RV):\n\n")
print(mt)

sink()

pdf.flnm <- "FCReadsFilter-Barplot.pdf"
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=10)
par(mfrow=c(2,1))
barplot(df$i.Reads,names.arg=1:nrow(df),las=2)
barplot(df$f.Reads,add=TRUE,col="lavender",yaxt="n")
title(main="Quality filter effect on raw reads")

par(mar=c(5,4.1,2,2)+0.1)
M <- t(data.matrix(res[,2:5]))
ymx <- max(apply(M,2,sum,na.rm=TRUE))*1.2
barplot(M,col=cls,ylim=c(0,ymx),xlab="file index",las=2)
legend("top",horiz=TRUE,fill=cls,legend=rownames(M),cex=0.6)

dev.off()

file.copy( file.path(repDir,txt.flnm), file.path(exportDir,txt.flnm),
           overwrite=TRUE )
file.copy( file.path(repDir,pdf.flnm), file.path(exportDir,pdf.flnm),
           overwrite=TRUE )
