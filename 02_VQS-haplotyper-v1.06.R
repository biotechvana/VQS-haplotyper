
###############################################################
###     VQS-HAPLOTYPES - QUASISPECIES DETECTON PIPELINE     ###
###############################################################

suppressMessages(library(optparse))

proj.nm <- "VQS_detection"

source("default_pars.R")

option_list <- list(
	make_option("--pmm_mx", action = "store", default = 3, type = "numeric", help = "Maximum number of mismatches in the specific primer [default %default]"),
	make_option("--min_reads", action = "store", default = 1, type = "numeric", help = "Minimum number of reads by sequences after repair [default %default]"),
	make_option("--max_ns", action = "store", default = 2, type = "numeric", help = "Maximum number of admissible Ns [default %default]"),
	make_option("--max_diffs", action = "store", default = 99, type = "numeric", help = "Maximum number of tolerated differences between a sequence and its reference [default %default]"),
	make_option("--max_gaps", action = "store", default = 50, type = "numeric", help = "Maximum number of admissible gaps [default %default]"),
	make_option("--ref_type", action = "store", default = "generic", type = "character", help = "Reference type to filter reads. Choose between generic or consensus [default %default]"),
	make_option("--method", action = "store", default = "Sum", type = "character", help = "Sum or Intersect. Sum takes the sum of the common haplotypes reads as a distribution, while Intersect takes the intersection as a distribution. [Default %default]"),
	make_option("--var_thr", action = "store", default = 0.5, type = "numeric", help = "Accept variant with an abundance above this value [default %default]"),
	make_option("--ab_thr", action = "store", default = 0.1, type = "numeric", help = "Second abundance filter for haplotypes [default %default]"),
	make_option("--min_size", action = "store", default = 1, type = "numeric", help = "Value between 0 and 1, which multiplies sequence length, to select minimum length of sequences after trimming adapters [default %default]"))

args <- parse_args(OptionParser(option_list = option_list))

pmm.mx <- args$pmm_mx
max.prdif <- pmm.mx
min.reads <- args$min_reads
max.Ns <- args$max_ns
max.diffs <- args$max_diffs
max.gaps <- args$max_gaps
ref.type <- args$ref_type
method <- args$method
var.thr <- args$var_thr
ab.thr <- args$ab_thr
min.size <- args$min_size

if ((method != "Sum") && (method != "Intersect")) {
	stop("Not recognized selected Method. Please, use Sum or Intersect.")
}

if ((ref.type != "generic") && (ref.type != "consensus")) {
	stop("Not recognized selected reference type. Please, use generic or consensus.")
}

tm <- integer(10)

cat("\nTrimming adaptors and primers\n")
print( (tm[1] = system.time( source("./R/RAV_FF_PrimersSplit_pl-v1.7.R") )[3]) ) 
cat("\nAligning and correcting haplotypes\n")
print( (tm[2] = system.time( source("./R/SeqFreqTable-MPFC-v4.5.R") )[3]) )
cat("\nFW and RV haplotype intersection\n")
print( (tm[3] = system.time( source("./R/nt-IntersectFW&RVHaplos-v5.3.R") )[3]) )
cat("\nAbundance filtering 0.5%\n")
print( (tm[4] = system.time( source("./R/AbFilterConsHaplos-050-v2.2.R") )[3]) )
cat("\nAbundance filtering 0.1%\n")
print( (tm[5] = system.time( source("./R/MinRdConsHaplos-010-v2.2.R") )[3]) )
cat("\nRare haplotypes profile\n")
print( (tm[6] = system.time( source("./R/RareHplProfile-v2.3.R") )[3]) )
cat("\nFinal coverage rport\n")
print( (tm[7] = system.time( source("./R/NtFastasSummary-v1.2.R") )[3]) )
print( (tm[8] = system.time( source("./R/Nt01FastasSummary-v1.2.R") )[3]) )

###  End of pipeline
cat("\nEnd of MiSeq nucleotide pipeline")
cat("\nElapsed globally: ",sum(tm),"\" (",round(sum(tm)/60,2),"\', ",
    round(sum(tm)/3600,2),"h)\n",sep="")

#---------------------------------------------------------------------------#
