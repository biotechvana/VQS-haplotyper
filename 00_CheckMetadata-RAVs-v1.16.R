
####################################################
###     Comprobació de fitxers i descriptors     ###
####################################################

library(Biostrings)
library(ShortRead)
library(data.table)

dataDir <- "./data"
runDir <- "./run"

cat("\nVerifying primer descriptors and files for pools in run")
cat("\n-------------------------------------------------------\n")

###  Llegim l'estructura de descripció de mostres
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Llista de pools en samples
pools <- unique(samples$Pool.Nm)

###  Llegim descriptors de primers
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Comprobar que tots els primers estan descrits
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
samples$Primer.idx <- idx
if( any(is.na(idx)) )
{  cat("No primers descriptor found for samples: ",
        paste(which(is.na(idx)),collapse=", "),"\n")
} else cat("\nAll primers in run are described.\n")  


###  Comprobar que hi hagi les RefSeqs de nt
RefSeqs <- as.character(
              readDNAStringSet(file.path(dataDir,"AmpliconRefSeqs.fna")))
rspatt <- paste(samples$Sbtp,".*_p",
                primers$FW.tpos[samples$Primer.idx],".",
                primers$RV.tpos[samples$Primer.idx],"$",sep="")
idx <- sapply(rspatt,function(patt) grep(patt,names(RefSeqs))[1])
if( any(is.na(idx)) )
{ cat("\nNo nucleotide reference sequence found for samples: ",
        paste(which(is.na(idx)),collapse=", "),"\n")
} else cat("\nAll nucleotide reference sequences found.\n")

###  Comprobar que hi hagi les RefSeqs d'aa
aaRefSeqs <- as.character(
              readAAStringSet(file.path(dataDir,"AmpliconAaRefSeqs.fna")))
rspatt <- paste("HAV.",
                primers$FW.tpos[samples$Primer.idx],".",
                primers$RV.tpos[samples$Primer.idx],"$",sep="")
idx <- sapply(rspatt,function(patt) grep(patt,names(aaRefSeqs))[1])
if( any(is.na(idx)) )
{ cat("\nNo amino acid reference sequence found for samples: ",
        paste(which(is.na(idx)),collapse=", "),"\n")
} else cat("\nAll amino acid reference sequences found.\n")
		
###  Llista de fastq en directori raw
flnms <- list.files(runDir)
if(length(flnms)==0)
  cat("\nNo fastq files in run folder.\n")
  
if(length(flnms))
{ snms <- sub("\\.fastq\\.gz$","",flnms)
  parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
  parts <- parts[,-c(3,5)]
  colnames(parts) <- c("PoolID","SmplID","Read")

  ###  Verificar que hi hagi un fastq R1 i R2 per cada pool en samples
  R1.IDs <- parts[parts[,3]=="R1","PoolID"]
  if( ! all(pools %in% R1.IDs) )
  { idx <- which(!(pools %in% R1.IDs)) 
    cat("\nNo R1 fastq file for pools:\n")
    print(pools[idx])
  } else cat("\nAll R1 fastq files in folder.\n")
  R2.IDs <- parts[parts[,3]=="R2","PoolID"]
  if( ! all(pools %in% R2.IDs) )
  { idx <- which(!(pools %in% R2.IDs)) 
    cat("\nNo R2 fastq file for pools:\n")
    print(pools[idx])
  } else cat("\nAll R2 fastq files in folder.\n")
}
cat("\n")
