
####################################################
###   Comprobación de ficheros y descripciones   ###
####################################################

library(Biostrings)
library(ShortRead)
library(data.table)

dataDir <- "./data"
runDir <- "./run"

cat("\nVerifying primer descriptors and files for pools in run")
cat("\n-------------------------------------------------------\n")

###  Lectura de la estructura de la descripción de las muestras
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Lista de pools en las muestras
pools <- unique(samples$Pool.Nm)

###  Lecturas de los descriptones de los adaptadores
primers <- fread(file.path(dataDir,"primers.csv"), sep="auto", header=T,
                      stringsAsFactors=F)
###  Comprobar que todos los primers estan descritos
idx <- sapply(samples$Primer.ID,function(pr) which(primers$Ampl.Nm==pr)[1])
samples$Primer.idx <- idx
if( any(is.na(idx)) )
{  cat("No primers descriptor found for samples: ",
        paste(which(is.na(idx)),collapse=", "),"\n")
} else cat("\nAll primers in run are described.\n")  


###  Comprobar la presencia de las secuencias de referencia (RefSeqs) en nucleótidos
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
		
###  Lista de fastq en directorio raw
flnms <- list.files(runDir)
if(length(flnms)==0)
  cat("\nNo fastq files in run folder.\n")
  
if(length(flnms))
{ snms <- sub("\\.fastq\\.gz$","",flnms)
  parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
  parts <- parts[,-c(3,5)]
  colnames(parts) <- c("PoolID","SmplID","Read")

  ###  Verificar que haya un fastq R1 y R2 por cada pool en samples
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
