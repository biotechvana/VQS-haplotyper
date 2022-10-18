##############################################################################
###                 MÓDULO: SeqFreqTable-FullCorrect          
###
###  Parte de los fasta de reads obtenidos de MiSeq una vez demultiplexados 
###    y recortados
###
###  Los haplotipos con más de max.gaps o max.Ns se eliminan
###  Si max.gaps es 0 se eliminan todas las reads con gaps
###  Si max.Ns es 0 se eliminan todos las reads con Ns
###
###  Se genera un informe sobre los resultados obtenidos.
##############################################################################

source("./R/global.v4.6.R")
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


###  Funciones de diferencias entre secuencias alineadas
#########################################################

##  Diferencias posición a posición
n.difs <- function(seq1,seq2)
{ sum( strsplit(seq1,split="")[[1]] != strsplit(seq2,split="")[[1]] ) }

##  Gaps a uno y otro lado
n.gaps <- function(seq1,seq2)
{ seq1 <- strsplit(seq1,split="")[[1]]
seq2 <- strsplit(seq2,split="")[[1]]
sum( (seq1=="-" & seq2!="-") | (seq1=="-" & seq2!="-") )
}

##  Diferencias sin tener en cuenta los gaps
n.edit <- function(seq1,seq2)
{ seq1 <- strsplit(seq1,split="")[[1]]
seq2 <- strsplit(seq2,split="")[[1]]
sum(seq1!=seq2)-sum((seq1=="-" & seq2!="-")|(seq1=="-" & seq2!="-"))
}

##  Porcentaje d eidntidad entre dos secuencias
pid2 <- function(seq1,seq2)
  100*(1-n.edit(seq1,seq2)/nchar(seq1))



###  Cuenta las Ns de cada secuencia
#####################################
has.Ns <- function(seq)
  as.vector(apply(as.matrix(seq),1,function(ch) sum(ch=="N")))


###  Determina las coincidencias de un segmento dentro de un conjunto de
###    fragmentos de secuencias.
#################################################################
procBySubstr <- function(seqs,patt,mm.mx,ambig)
{
  ##  Buscar coincidencias.
  matches <- vmatchPattern(patt,subject=seqs,fixed=!ambig,max.mismatch=mm.mx,
                           with.indels=TRUE) ## CAMBIO. with.indels = TRUE
  ##  elementLengths() es de dimensión completa con las posiciones de
  ##    no coincidencia a 0. Evitar coincidencias múltiples.
  flags <- elementLengths(matches)==1
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


###  Corregir las Ns de secuencias con menos de max.Ns
###   Si hay más de max.Ns en una secuencia la deja inalterada.
####################################################################
CorrNs <- function(cseqs, RefSeq)
{
  seqlist <- list()
  nslist <- list()
  for (i in 1:length(as.matrix(cseqs))) {
    ## Secuencias a corregir (hasta un máximo de max.Ns)
    ns <- has.Ns(DNAStringSet(cseqs[i]))
    flc <- ns>0 & ns<=max.Ns
    if (sum(flc)==0) {
      seqlist[[i]] <- as.character(cseqs[i])
      nslist[[i]] <- as.integer(ns)
    } else {
      ## Posición a corregir
      sqs <- DNAStringSet(cseqs[i][flc])
      cm <- consensusMatrix(sqs)
      pos <- which(cm["N",]>0)
      
      ## Base en esta posición por la RefSeq
      let <- strsplit(RefSeq[i], split="")[[1]][pos]
      for (j in 1:length(pos)) {
        ## vector de lógicos indicando las secuencias que tienen una N en esta posición
        fl <- letter(as.character(sqs), pos[j]) == "N"
        ## matriz de lógicos indicando las posiciones donde hay Ns a corregir
        at <- matrix(FALSE, nrow=sum(fl), ncol=nchar(sqs[1]))
        at[cbind(1:nrow(at), pos[j])] <- TRUE
        ## corrección
        sqs[fl] <- replaceLetterAt(DNAStringSet(sqs[fl]), at=at, letter=rep(let[j], sum(fl)))
      }
      seqlist[[i]] <- as.character(sqs)
      nslist[[i]] <- as.integer(ns)
    }
  }
  list(cseqs=as.character(seqlist),ns=as.integer(nslist))
}


###  Da nombre a los haplotipos en función del número de mutaciones
###  respecto la referencia (que puede ser la máster) y de su frecuencia
###  poblacional, y guarda las secuencias en un fichero fasta.
########################################################################
SaveHaplotypes <- function(ifl,bseqs,nr)
{
  ##  Nombre genérico del haplotipo
  code <- paste(FlTbl$Pat.ID[ifl],FlTbl$Ampl.Nm[ifl],FlTbl$Str[ifl],sep=".")
  
  ##  Determinar diferencias respecto a la máster
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[1])
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
  nms <- paste(code,nm,zeroFillInt2Char(isq,4),sep=".")
  
  ##  Cabecera fasta con nombre, núemro de reads y frecuencia relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")
  
  ##  Guardar a fasta
  flnm <- sub("TCA.454Reads","nt.FCHaplos",FlTbl$File.Name[ifl])
  write.fasta(DNAStringSet(bseqs),file.path(filtDir,flnm))
  
  list(bseqs=bseqs,nr=nr,nm=nm)
}


###  Determinar máster.
########################
CorrectMaster <- function(bseqs,nr,RefSeq=NULL)
{
  ##  Si no hay reads
  if(!length(bseqs))
  { cat("\nNo reprentative sequences found.\n")
    cat("\n-------------------------------------------------\n\n")
    return(NULL)  ## Conjunto vacío.
  }
  
  ##  Si no hay secuencia de referencia, devolver
  if(is.null(RefSeq)) 
    return( list(bseqs=bseqs, nr=nr, mseq=bseqs[1]) )

  list(bseqs=bseqs, nr=nr, mseq=bseqs[1])
}


###  Alinea haplotipos respecto a la máster o referencia y los filtra
###   y corrige según los parámetros del programa
####################################################################
FilterCorrectHaplos <- function(ifl,bseqs,nr,mnrd=2,RefSeq=NULL)
{
  low.freq <- 0
  res <- data.frame(tmdifs=integer(1),tmNs=integer(1),tmdels=integer(1),
                    corr=integer(1),f.reads=integer(1),f.haplo=integer(1))

  cons.seq <- ""
  names(cons.seq) <- paste(FlTbl$Pat.ID[ifl],FlTbl$Ampl.Nm[ifl],FlTbl$Str[ifl],
                           sep=".")
  
  ##  Alinear respecto a la RefSeq o la máster
  mseq <- RefSeq
  if(is.null(mseq)) mseq <- bseqs[1]  ##  usar la master
  psa <- pairwiseAlignment(pattern=bseqs,subject=mseq)
  
  
  ##  Invocar aligned() es fundamental: devuelve las secuencias con gaps donde
  ##   hay inserciones respecto a la referencia, y con las posiciones que dan
  ##   gaps en la ref eliminadas, lo que comportará gaps a final de cadena.
  bseqs <- as.character(pattern(psa))
  refs_aln <- as.character(subject(psa))
  
  ##  Informar de secuencias con demasiadas diferencias
  ndiffs <- nmismatch(psa)
  dif.fl <- ndiffs <= max.diffs         
  htmd <- sum(!dif.fl)
  res$tmdifs <- rtmd <- sum(nr[!dif.fl])    ###
  if(htmd)
  { cat("\nExcluded haplotypes with too many differences:",sum(!dif.fl))
    cat("  (",rtmd," reads)",sep="")
  }
  
  ##  Informar de deleciones (gaps en las reads)
  dels <- deletion(psa)
  del.fl <- sapply(dels,function(l) length(l)>0)
  if(sum(del.fl))
  { cat("\nHaplotypes with deletions:",sum(del.fl) )
    cat("  (",sum(nr[del.fl])," reads)",sep="")
  }
  
  ##  Informar de inserciones (gaps en la referencia)
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
  { cNs <- sum(fl.Ns.corr) # número de haplotipos corregidos
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
  
  ##  Detectar número de gaps
  fl.tmg <- rep(FALSE,length(bseqs))
  if(sum(del.fl)+sum(ins.fl))
  { 
    ngap <- alphabetFrequency(DNAStringSet(bseqs))[,"-"]
    fl.ng <- ngap > 0 & ngap <= max.gaps
    fl.tmg <- ngap > max.gaps
  }
  
  cat("\n\nTotal haplotypes excluded by Ns, gaps or differences:",
      sum( fl.tmn | fl.tmg | !dif.fl ) )
  cat("  (",sum(nr[fl.tmn | fl.tmg | !dif.fl])," reads)",sep="")
  
  ##  Eliminar los haplotipos que no pasan los filtros
  bseqs <- bseqs[!fl.tmn & !fl.tmg & dif.fl]
  nr <- nr[!fl.tmn & !fl.tmg & dif.fl]
  
  ##  Si no hay secuencias válidas
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
    
    o <- order(nr,decreasing=TRUE)
    bseqs <- bseqs[o]
    nr <- nr[o]
    if( mseq != bseqs[1] )
    { cat("\n  Warning: a new master emerged after the repairs")
      cat("\n           with frequency",nr[1],"reads")
      mseq <- bseqs[1]
    }
    
    ##  Post-filtro por número de reads
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
      return(list(res=res,mseq=NULL,nrmseq=0,low.freq=low.freq))  ## Vacío.
    }
    
    ##  Generar la consenso poblacional.
    consMat <- SeqsTbl.w(bseqs,nr)
    cons.seq <- ConsSeq(consMat)
    names(cons.seq) <- "000|0"
    
    ##  Diferencias entre la máster y la consenso 
    cmdifs <- n.edit(bseqs[1],cons.seq)
    cat("\n\nConsensus – master differences:",cmdifs)
  }
  cat("\nFinal reads for the master sequence:",nr[1])
  
  ## Ordenar por mutaciones y frecuencias, poner nombre y guardar haplotipos
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

###  Lectura de la estructura de las descripciones de las muestras
##################################################
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)

###  Localizar amplicones del experimento
##################################################
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
if( any(is.na(idx)) )
  stop("No primers descriptor found for samples: ",
       paste(which(is.na(idx)),collapse=", "))
samples$Primer.idx <- idx

###  Cargar secuencias de referencia genéricas por cada muestra
#################################################################
RefSeqs <- as.character(
  readDNAStringSet(file.path(dataDir,"AmpliconRefSeqs.fna")))
idx <- sapply(samples$RefSeq.ID,function(pr) which(names(RefSeqs)==pr)[1])
if( any(is.na(idx)) )
  stop("No reference sequence found for samples: ",
       paste(which(is.na(idx)),collapse=", "))
###  Vincular una RefSeq a cada muestra
samples$RefSeq.idx <- idx

##  Lectura de la tabla de ficheros a tratar
#############################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
FlTbl$RefSeq.idx <- 0
for(i in 1:nrow(FlTbl))
{ k <- which(samples$Patient.ID == FlTbl$Pat.ID[i] &
               samples$Primer.ID == FlTbl$Ampl.Nm[i])
FlTbl$RefSeq.idx[i] <- samples$RefSeq.idx[k]
}  

##  Inicio del informe
#######################
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

res <- data.frame(i.reads=integer(),low.freq=integer(),
                  tmdifs=integer(),tmNs=integer(),tmdels=integer(),
                  corr=integer(),f.reads=integer(),f.haplo=integer())

##  Loop sobre ficheros a analitzar
#######################################
ampl.set <- unique(sort(FlTbl$Ampl.Nm))
mseqs <- character(nrow(FlTbl))
for( i in 1:nrow(FlTbl) )
{
  ##  Secuencia de referencia para las reads limpias del fichero
  RefSeq <- NULL
  pr.id <- FlTbl$Pr.ID[i]
  if( !is.null(RefSeqs) ) 
    RefSeq <- RefSeqs[FlTbl$RefSeq.idx[i]]
  
  ##  Leer fichero fasta
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
  
  ##  Identificar máster
  rseq <- NULL
  if(!is.null(RefSeq) & ref.type=="generic") rseq <- RefSeq
  lst <- CorrectMaster(bseqs,nr,rseq)
  if(is.null(lst)) 
    next
  ##  Guardarla
  mseqs[i] <- lst$mseq
  
  ##  Filtrar y corregir haplotipos por la máster si no hay ninguna RefSeq, 
  ##    por la máster si la RefSeq es "generic", o
  ##    por la misma RefSeq si es "consensus"
  rseq <- NULL
  if(ref.type=="consensus") rseq <- RefSeq
  lst <- FilterCorrectHaplos(i,lst$bseqs,lst$nr,mnrd=min.reads,RefSeq=rseq)
  res[i,"low.freq"] <- lst$low.freq
  res[i,3:8] <- lst$res
  if(lst$res$f.reads > 0 & !is.null(lst$mseq))
    mseqs[i] <- lst$mseq
  
}      ##  Final de loop sobre amplicones

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

##  Muestras perdidas.
nlst <- sum(df$i.Reads>0 & df$f.Reads==0)
if(nlst)
{ cat("\n\nSamples lost in the filters:  ",nlst," (",
      round(nlst/nrow(df)*100,2),"%)\n",sep="")
  fl <- df$f.Reads>0
  cat("\nGlobal table by files:\n")
  print(table(df$Ampl.Nm[fl],df$Pat.ID[fl]))
}

##  Tabla de reads por amplicón (FW+RV)
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
