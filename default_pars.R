
###############################################
###  PARAMETRES DEL PIPE-LINE DE nt en HCV
###############################################

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
resultsDir <- resDir <- "./results"
exportDir <- expDir <- "./export"
tempDir <- "./tmp"
ntDir <- "./nt"
nt01Dir <- "./nt.01"
aaDir <- "./aa"

#-----------------------------------------------------------------------#

###  Per splitByPrimers  ###

##  TRIM PRIMERS
## Definim parametres per separació de reads:
pmm.mx <- 3        ## Nombre màxim de mismatch en el primer específic
max.prdif <- pmm.mx
min.len <- 200     ## Longitud mínima per considerar una seqüència
###  Depèn que hi hagi adaptadors 454, MIDs i/o M13
target.io <- 1   # 10   25   50
target.in <- 30  # 55   55   80

#-----------------------------------------------------------------------#

## Definim paràmetres de filtrat de reads
min.reads <- 1         ## Nombre mínim de reads per seqüència, post reparacions.
max.Ns <- 2            ## Nombre màxim de Ns admisibles.
max.diffs <- 99        ## Nombre màxim de diferències tolerat.
max.gaps <- 3          ## Nombre màxim de gaps admisibles.
ref.type <- "generic"  ## Un d'entre "generic" o "consensus"

#-----------------------------------------------------------------------#

## Definim paràmetres per la intersecció d'haplotips

###  Mètode de confrontació FW i RV
###  'Sum'  pren la suma de reads d'haplotips comuns com a distribució
###  'Intersect' pren la intersecció com a distribució
method <- "Sum"    ###  Un de c("Intersect","Sum")
min.rd <-   1      ###  Si cal filtrat per mínim nombre de reads
a.cut  <-  0.00    ###  Llindar d'abundància previ a la intersecció en %
ni.thr <-   2      ###  Llindar d'abundància per salvar no comuns en %

#-----------------------------------------------------------------------#

## Definim paràmetres per acceptar variants

var.thr <- 0.50   ###  Acceptar variants per sobre del 0.5%

#-----------------------------------------------------------------------#
