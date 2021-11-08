
##############################################
###     MiSeq - HCV SUBTYPING PIPELINE     ###
###      Quality assessment pipeline       ###
##############################################

library(Biostrings)
library(ShortRead)

proj.nm <- "VHIR38_CBM"

codeDir <- "./R"
dataDir <- "./data"
runDir <- "./run"
flashDir <- "./flash"
flashFiltDir <- "./flashFilt"
repDir <- "./reports"

min.len <- 200
min.ov <- 20     # Minimum overlap between R1 and R2
max.ov <- 300    # Maximum overlap between R1 and R2
err.lv <- 0.10   # Fraction of accepted mismatches in overlap 
flash <- "./FLASH-1.2.11/flash"
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

###  Hi ha M13 + primer específic, hi poden haver MIDs també
target.io <- 20-5
target.in <- target.io+50

tm <- integer(10)

#cat("\nRunning quality assessment on R1 and R2 fastq files\n")
#print( (tm[1] = system.time( source("./R/R1R2_LoqQ2N_pl-v1.27.R") )[3]) ) 
cat("\nRunning FLASH to extend reads\n")
print( (tm[2] = system.time( source("./R/R1R2_to_FLASH_pl.R") )[3]) ) 
cat("\nRunning QC by position, R1, R2 and FLASH fastq files\n")
print( (tm[3] = system.time( source("./R/PoolQCbyPos-v2.3.R") )[3]) ) 
cat("\nRunning QC by read on FLASH fastq files\n")
print( (tm[4] = system.time( source("./R/PoolQCbyRead-v1.2.R") )[3]) ) 
cat("\nRunning length peaks on FLASH fastq files\n")
print( (tm[5] = system.time( source("./R/LenPeaksByPool-v2.R") )[3]) ) 

##  Accept max of 5% bases below Q30 by read
ThrQ30 <- 0.05
cat("\nFiltering haplotypes by Q30\n")
print( (tm[6] = system.time( source("./R/FiltByQ30-v1.1.R") )[3]) ) 
cat("\nRunning QC by position after filtering fastq files\n")
print( (tm[7] = system.time( source("./R/PoolFiltQCbyPos-v2.2.R") )[3]) ) 

###  End of pipeline
cat("\nEnd of MiSeq quality assessment pipeline")
cat("\nElapsed globally: ",sum(tm),"\" (",round(sum(tm)/60,2),"\', ",
    round(sum(tm)/3600,2),"h)\n",sep="")
	
#---------------------------------------------------------------------------#
	
tm
 # [1]   0.00 533.77 708.37 413.05  98.87
 
# End of MiSeq quality assessment pipeline
# Elapsed globally: 1754.06" (29.23', 0.49h)
