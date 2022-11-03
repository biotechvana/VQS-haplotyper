
################################################
###   PARÁMETROS DEL PIPELINE VQS-HAPLOTYPER   #
################################################

stopifnot(suppressMessages(require(Biostrings)))
stopifnot(suppressMessages(require(stringr)))
stopifnot(suppressMessages(require(RColorBrewer)))

dataDir <- "./data"
codeDir <- "./R"
runDir <- "./run"
flashDir <- "./flash"
flashFiltDir <- "./flashFilt"
splitDir <- "./splits"
trimDir <- "./trim"
filtDir <- "./filt"
joinDir <- "./join"
repDir <- "./reports"
resultsDir <- "./results"
exportDir <- "./export"
tempDir <- "./tmp"
ntDir <- "./nt"
nt01Dir <- "./nt.01"
aaDir <- "./aa"

#-----------------------------------------------------------------------#

###  Para splitByPrimers  ###

##  TRIM PRIMERS
###  Depende de si hay adaptadores 454, MIDs y/o M13
target.io <- 1   # 10   25   50
target.in <- 30  # 55   55   80

#-----------------------------------------------------------------------#
