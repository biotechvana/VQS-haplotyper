
##############################################
###      MiSeq -  SUBTYPING PIPELINE       ###
###      Quality assessment pipeline       ###
##############################################

library(Biostrings)
library(ShortRead)

suppressMessages(library(optparse))

option_list <- list(
	make_option("--min_len", action = "store", default = 200, type = "numeric", help = "Minimum length to consider a sequence [default %default]"),
	make_option("--min_ov", action = "store", default = 20, type = "numeric", help = "Minimum overlap between R1 and R2 in FLASH [default %default]"),
	make_option("--max_ov", action = "store", default = 300, type = "numeric", help = "Maximum overlap between R1 and R2 in FLASH [default %dfault]"),
	make_option("--err_lv", action = "store", default = 0.10, type = "numeric", help = "Fraction of accepted mismatches in overlap in FLASH [default %default]"),
	make_option("--thrQ30", action = "store", default = 0.05, type = "numeric", help = "Maximum percentage of bases below Q30 accepted by read [default %default]"))

proj.nm <- "VQS_detection"

codeDir <- "./R"
dataDir <- "./data"
runDir <- "./run"
flashDir <- "./flash"
flashFiltDir <- "./flashFilt"
repDir <- "./reports"

min.len <- args$min_len   # Minimum length to consider a sequence
min.ov <- args$min_ov   # Minimum overlap between R1 and R2
max.ov <- args$max_ov   # Maximum overlap between R1 and R2
err.lv <- args$err_lv   # Fraction of accepted mismatches in overlap


if(.Platform$OS.type == "unix") {
	flash <- "./FLASH-1.2.11/flash"
	} else {
	flash <- "./FLASH-1.2.11/windows/flash.exe"
	}
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

###  Hay M13 + primer específico, puede haber MIDs también
target.io <- 20-5
target.in <- target.io+50

tm <- integer(10)

cat("\nRunning FLASH to extend reads\n")
print( (tm[2] = system.time( source("./R/R1R2_to_FLASH_pl.R") )[3]) ) 
cat("\nRunning QC by position, R1, R2 and FLASH fastq files\n")
print( (tm[3] = system.time( source("./R/PoolQCbyPos-v2.3.R") )[3]) ) 
cat("\nRunning QC by read on FLASH fastq files\n")
print( (tm[4] = system.time( source("./R/PoolQCbyRead-v1.2.R") )[3]) ) 
cat("\nRunning length peaks on FLASH fastq files\n")
print( (tm[5] = system.time( source("./R/LenPeaksByPool-v2.R") )[3]) ) 

##  Accept max of 5% bases below Q30 by read
ThrQ30 <- args$thrQ30
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
