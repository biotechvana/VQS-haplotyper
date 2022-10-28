
dataDir <- "./data"
repDir <- "./reports"
library(data.table)

###  Cargar variables rare.w, Tnr, Tnh, log.cuts
load(file.path(repDir,"RareHplWeight.RData"))
###  Cargar descriptores de las muestras
samples <- fread(file.path(dataDir,"samples.csv"), sep="auto", header=T,
                 stringsAsFactors=F)
### Sincronizar muestras con las otras estructuras
o <- sapply(rownames(rare.w),function(nm) 
  which(paste(samples$Patient.ID,samples$Primer.ID,sep=".")==nm))
samples <- samples[o,]

###  Factores
ampl <- unique(samples$Primer.ID)
trt <- unique(samples$Treat)
pbase <- unique(samples$pBase)
ptrt <- unique(samples$pTreat)

pdf.flnm <- file.path(repDir,"RareHaplotypesProfile.pdf")
pdf(pdf.flnm,paper="a4",width=6.5,height=10.5)
pal <- c("navy","maroon","darkgreen")

par(mfrow=c(3,1))

txt.flnm <- file.path(repDir,"RareHaplotypesWeight.txt")
sink(txt.flnm)

lims <- log10(c(0.0001,0.001,0.01))
jdx <- sapply(lims,function(lm) which(log.cuts==lm))
df.rh <- data.frame(AllHpl=Tnh,AllReads=Tnr)
m <- matrix(nrow = length(unique(samples$Primer.ID)), ncol = 2)

for (i in 1:length(unique(samples$Primer.ID))) {
  idx <- which(samples$Primer.ID==unique(samples$Primer.ID[i]))
  M <- rare.w[idx,jdx]
  seM <- sqrt(M*(1-M)/Tnr[idx])
  M <- M*100
  seM <- seM*100
  if (length(idx) == 1) {
    M <- as.data.frame(matrix(M,ncol = length(M),byrow = T))
    seM <- as.data.frame(matrix(seM,ncol = length(seM),byrow = T))
    colnames(M) <- colnames(seM) <- c("0.01%","0.1%","1%")
    rownames(M) <- rownames(seM) <- rownames(rare.w)[[idx]]
  } else{
    colnames(M) <- colnames(seM) <- c("0.01%","0.1%","1%")
  }

  if (length(idx) == length(unique(samples$Patient.ID))) {
    longest_idx <- idx
    longest_seM <- seM
  }
  
  if (i == "1") {
    o <- order(M[,"1%"])
  }
  
  nam <- paste("M", i, sep = "")
  M_assign <- M
  rownames(M_assign) <- gsub("\\..*","",rownames(M_assign))
  assign(nam, M_assign)
  m[i,1] <- nam
  m[i,2] <- length(get(nam)[,1])
  nms <- paste(samples$Patient.ID[idx],samples$Primer.ID[idx],sep=".")
  
  ylm <- c(0,70)
  plot(1:nrow(M),M[,1],type="l",font=2,col=pal[1],xaxt="n",xlab="",
       ylab="Weight (%)",ylim=ylm)
  pit <- sapply(2:ncol(M),function(j) 
    lines(1:nrow(M),M[,j],col=pal[j]))
  pit <- sapply(1:ncol(M),function(j) 
    arrows(x0=1:nrow(M),y0=M[,j]-2*seM[,j],y1=M[,j]+2*seM[,j],col=pal[j],
           code=3,angle=90,len=0.04))
  abline(v=1:nrow(M),lty=4,col="gray")
  axis(1,at=1:nrow(M),lab=nms,cex.axis=0.8,las=2)
  title_plot <- paste("Rare haplotypes weight -", unique(samples$Primer.ID)[i], sep = " ")
  title(main=title_plot)
  legend("top",horiz=TRUE,lwd=3,col=pal,bg="white",legend=colnames(M),cex=0.8)
  
  cat(sprintf("\nAmplicon %s:\n", unique(samples$Primer.ID)[i]))
  cat("\nAll reads and haplotypes\n\n")
  print(df.rh[idx,])
  cat("\nRare haplotypes weight (%)\n\n")
  print(round(M,2))
  cat("\nEstandard errors\n\n")
  print(round(seM,4))
  
}

m <- m[order(m[,2], decreasing = TRUE),]

###  Mean of the amplicons
addition = NULL
for (i in 1:length(m[,1])) {
  if (i == 1) {
    addition = get(m[i])
  } else {
    merged <- merge(addition, get(m[i]), by = 0, all = TRUE)
    rownames(merged) <- merged[,1]
    merged <- merged[,2:ncol(merged)]
    merged[,1:ncol(merged)] = apply(merged[,1:ncol(merged)], 2, function(x) as.numeric(as.character(x)))
    colnames(merged) <- gsub("\\.[x|y]","",colnames(merged))
    merged[is.na(merged)] <- 0
    addition <- sapply(unique(colnames(merged)), function(x) rowSums(merged[, colnames(merged)==x, drop=FALSE]))
  }
}

M <- addition/length(m[,1])
plot(1:nrow(M),M[,1],type="l",font=2,col=pal[1],xaxt="n",xlab="",
     ylab="Weight (%)",ylim=ylm)
pit <- sapply(2:ncol(M),function(j) 
  lines(1:nrow(M),M[,j],col=pal[j]))
pit <- sapply(1:ncol(M),function(j) 
  arrows(x0=1:nrow(M),y0=M[,j]-longest_seM[,j],y1=M[,j]+longest_seM[,j],col=pal[j],
         code=3,angle=90,len=0.02))
abline(v=1:nrow(M),lty=4,col="gray")
nms <- paste(samples$pBase[longest_idx],samples$Treat[longest_idx],samples$pTreat[longest_idx],sep=".")
axis(1,at=1:nrow(M),lab=nms,cex.axis=0.8,las=2)
title(main="Rare haplotypes weight - Mean of amplicons")
legend("top",horiz=TRUE,lwd=3,col=pal,bg="white",legend=colnames(M),cex=0.8)

cat("\nMean of the amplicons:\n")
cat("\nRare haplotypes weight (%)\n\n")
print(round(M,2))

sink()     
dev.off()
